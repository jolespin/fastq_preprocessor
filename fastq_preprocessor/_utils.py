import os
import sh
from loguru import logger


def create_directory(path):
    os.makedirs(path, exist_ok=True)
    return path


def format_path(path):
    return os.path.realpath(os.path.expandvars(os.path.expanduser(path)))


def assert_acceptable_arguments(value, acceptable):
    assert value in acceptable, "{!r} is not one of the acceptable options: {}".format(value, acceptable)


class Pipeline:
    """Lightweight pipeline runner using sh and loguru, replacing genopype.ExecutablePipeline."""

    def __init__(self, name, description, f_cmds, checkpoint_directory, log_directory):
        self.name = name
        self.description = description
        self.f_cmds = f_cmds
        self.checkpoint_directory = checkpoint_directory
        self.log_directory = log_directory
        self.steps = []

    def add_step(self, id, description, step, cmd, input_filepaths, output_filepaths,
                 validate_inputs=True, validate_outputs=True):
        self.steps.append({
            "id": id,
            "description": description,
            "step": step,
            "cmd": cmd,
            "input_filepaths": input_filepaths,
            "output_filepaths": output_filepaths,
            "validate_inputs": validate_inputs,
            "validate_outputs": validate_outputs,
        })

    def _cmd_to_str(self, cmd):
        parts = []
        for token in cmd:
            token = str(token).strip()
            if token:
                parts.append(token)
        return " ".join(parts)

    def compile(self):
        """Write all step commands to the commands file."""
        for step_info in self.steps:
            cmd_str = self._cmd_to_str(step_info["cmd"])
            self.f_cmds.write("# Step {}: {} | {}\n".format(
                step_info["step"], step_info["id"], step_info["description"]))
            self.f_cmds.write(cmd_str + "\n\n")

    def execute(self, restart_from_checkpoint=None):
        """Execute all pipeline steps, optionally restarting from a checkpoint."""
        for step_info in self.steps:
            step_num = step_info["step"]
            step_id = step_info["id"]
            checkpoint_file = os.path.join(self.checkpoint_directory, "{}__{}" .format(step_num, step_id))

            # Skip steps already completed before the restart point
            if restart_from_checkpoint is not None and step_num < restart_from_checkpoint:
                if os.path.exists(checkpoint_file):
                    logger.info("[Step {}] Skipping ({}): checkpoint exists", step_num, step_id)
                    continue

            logger.info("[Step {}] {} ({})", step_num, step_info["description"], step_id)

            # Validate inputs
            if step_info["validate_inputs"]:
                for fp in step_info["input_filepaths"]:
                    if "*" not in fp and not os.path.exists(fp):
                        raise FileNotFoundError("[Step {}] Input not found: {}".format(step_num, fp))

            cmd_str = self._cmd_to_str(step_info["cmd"])
            logger.info("[Step {}] Command: {}", step_num, cmd_str)
            log_base = os.path.join(self.log_directory, "{}__{}" .format(step_num, step_id))

            returncode = 0
            with open("{}.stdout".format(log_base), "w") as stdout_fh, \
                 open("{}.stderr".format(log_base), "w") as stderr_fh:
                try:
                    sh.bash("-c", cmd_str, _out=stdout_fh, _err=stderr_fh)
                except sh.ErrorReturnCode as e:
                    returncode = e.exit_code
                    with open("{}.returncode".format(log_base), "w") as f:
                        f.write("{}\n".format(returncode))
                    raise RuntimeError(
                        "[Step {}] {} failed (exit code {}). Check: {}.stderr".format(
                            step_num, step_id, returncode, log_base)
                    )

            with open("{}.returncode".format(log_base), "w") as f:
                f.write("{}\n".format(returncode))

            # Validate outputs
            if step_info["validate_outputs"]:
                for fp in step_info["output_filepaths"]:
                    if "*" not in fp and not os.path.exists(fp):
                        raise FileNotFoundError("[Step {}] Output not found: {}".format(step_num, fp))

            # Mark checkpoint
            open(checkpoint_file, "w").close()
            logger.success("[Step {}] Completed ({})", step_num, step_id)
