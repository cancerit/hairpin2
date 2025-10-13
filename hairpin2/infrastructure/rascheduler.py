from collections.abc import Iterable, Mapping
from typing import Any, Self, cast

from pydantic import ValidationError

from hairpin2.infrastructure.process import ReadAwareProcess, ReadAwareProcessProtocol
from hairpin2.infrastructure.process_engines import ProcessKindEnum
from hairpin2.infrastructure.process_params import FixedParams, RunParams
from hairpin2.infrastructure.structures import FlagResult


class ConfigError(Exception):
    pass


class ExecutorError(RuntimeError):
    pass


class ProcessError(TypeError):
    pass


class RAExec:
    _mandate_excludes: bool
    _taggers: tuple[ReadAwareProcessProtocol[RunParams, FixedParams | None, None], ...]
    _flaggers: tuple[ReadAwareProcessProtocol[RunParams, FixedParams | None, FlagResult], ...]

    def __init__(
        self,
        taggers_in_order: Iterable[ReadAwareProcessProtocol[RunParams, FixedParams | None, None]],
        flaggers_in_order: Iterable[
            ReadAwareProcessProtocol[RunParams, FixedParams | None, FlagResult]
        ],
        mandate_excludes: bool,
    ):
        self._mandate_excludes = mandate_excludes
        _, tags_added = self.validate_tagger_order(
            taggers_in_order, mandate_excludes, raise_on_fail=True
        )
        _ = self.validate_flagger_exec(
            flaggers_in_order, tags_added, mandate_excludes, raise_on_fail=True
        )
        self._taggers = tuple(taggers_in_order)
        self._flaggers = tuple(flaggers_in_order)

    @property
    def mandate_excludes(self):
        return self._mandate_excludes

    @property
    def taggers(self):
        return self._taggers

    @property
    def flaggers(self):
        return self._flaggers

    @staticmethod
    def _partition_processes(
        processes: Iterable[ReadAwareProcessProtocol[Any, Any, Any]], raise_on_fail: bool = False
    ) -> tuple[
        bool,
        tuple[
            tuple[ReadAwareProcessProtocol[RunParams, FixedParams | None, None], ...],
            tuple[ReadAwareProcessProtocol[RunParams, FixedParams | None, FlagResult], ...],
        ],
    ]:
        taggerl: list[ReadAwareProcessProtocol[RunParams, FixedParams | None, None]] = []
        flaggerl: list[ReadAwareProcessProtocol[RunParams, FixedParams | None, FlagResult]] = []
        taggers_all_registered = False
        valid = True
        for proc in processes:
            match proc.ProcessType:
                # TODO: this should use a runtime check rather than a cast
                case ProcessKindEnum.TAGGER:
                    if taggers_all_registered:
                        if raise_on_fail:
                            raise ExecutorError("tagger appeared in execution order after flagger")
                        valid = False
                    taggerl.append(
                        cast(ReadAwareProcessProtocol[RunParams, FixedParams | None, None], proc)
                    )
                case ProcessKindEnum.FLAGGER:
                    if not taggers_all_registered:
                        taggers_all_registered = True
                    flaggerl.append(
                        cast(ReadAwareProcessProtocol[RunParams, FixedParams, FlagResult], proc)
                    )
                case _:
                    raise ProcessError(
                        f"Unrecognised process type found on process {proc.ProcessNamespace}"
                    )
        return valid, (tuple(taggerl), tuple(flaggerl))

    @staticmethod
    def check_namespacing_clashfree(processes: Iterable[ReadAwareProcessProtocol[Any, Any, Any]]):
        nsl = [proc.ProcessNamespace for proc in processes]
        return len(nsl) == len(set(nsl))  # clashes if false

    @staticmethod
    def validate_exec(
        processes: Iterable[ReadAwareProcessProtocol[Any, Any, Any]],
        mandate_excludes: bool,
        raise_on_fail: bool,
    ) -> tuple[
        bool,
        tuple[
            tuple[ReadAwareProcessProtocol[Any, Any, None], ...],
            tuple[ReadAwareProcessProtocol[Any, Any, FlagResult], ...],
        ],
    ]:
        valid, part_procs = RAExec._partition_processes(processes, raise_on_fail)
        if valid:
            valid, tags_added = RAExec.validate_tagger_order(
                part_procs[0], mandate_excludes, raise_on_fail
            )
            if valid:
                valid = RAExec.validate_flagger_exec(
                    part_procs[1], tags_added, mandate_excludes, raise_on_fail
                )
        return valid, part_procs

    # make wrapped read register itself as marked in its parent readview if there is one
    @staticmethod
    def validate_tagger_order(
        processes: Iterable[ReadAwareProcessProtocol[RunParams, FixedParams | None, None]],
        mandate_excludes: bool,
        raise_on_fail: bool = False,
    ) -> tuple[bool, set[str]]:
        """
        Validate that the order of tagger processes as provided in the processes iterable can be run without error.
        In practice this means checking that at a minimum, all marks/tags listed in a given process as required (i.e.
        in require_marks) have been checked and applied by processes occurring prior in the run order.

        Args:
            processes (Iterable): Iterable of processes, the order of which defines the run order
            mandate_excludes (bool): if true, then for each process, all marks listed as exclude_marks must have been applied by tagger processes coming prior in the run order, in addition to require_marks
            raise_on_fail (bool): whether to throw an error on receiving an illegal order, or just return a bool
        """
        tags_set: set[str] = set()
        valid = True
        for proc in processes:
            assert proc.AddsMarks is not None
            reqs_fulfilled = proc.require_marks.issubset(tags_set)
            if mandate_excludes:
                reqs_fulfilled = reqs_fulfilled and proc.exclude_marks.issubset(tags_set)
            if not reqs_fulfilled:
                if raise_on_fail:
                    raise ExecutorError("Invalid tagger execution order")
                valid = False
                break
            if tags_set & proc.AddsMarks:  # if overlap
                if raise_on_fail:
                    raise ExecutorError("taggers appear to add same tag. This is not supported")
                valid = False
                break
            tags_set = tags_set | proc.AddsMarks
        return (valid, tags_set) if valid else (valid, set())

    @staticmethod
    def validate_flagger_exec(
        processes: Iterable[ReadAwareProcessProtocol[RunParams, FixedParams | None, FlagResult]],
        tags_set: set[str],
        mandate_excludes: bool,
        raise_on_fail: bool = False,
    ):
        valid = True
        for proc in processes:
            reqs_fulfilled = proc.require_marks.issubset(tags_set)
            if mandate_excludes:
                reqs_fulfilled = reqs_fulfilled and proc.exclude_marks.issubset(tags_set)
            if not reqs_fulfilled:
                if raise_on_fail:
                    raise ExecutorError("Invalid tagger execution order")
                valid = False
                break
        return valid

    @classmethod
    def from_config(
        cls,
        configd: dict[str, Any],
        proc_pool: Iterable[
            type[ReadAwareProcessProtocol[RunParams, FixedParams | None, FlagResult | None]]
        ],
    ) -> Self:
        """
        Initialise read-aware executor from config dict.

        Args:
            configd (dict[str, Any]): Dictionary describing how processes should execute.
            proc_pool (Iterable[type[ReadAwareProcess]]): Iterable of process definitions which may be instatiated for a given run.
        """
        for proc in proc_pool:
            if not issubclass(proc, ReadAwareProcess):
                raise ProcessError(
                    f"Process {proc} is not a concrete implementation (subclass) of abstract base class ReadAwareProcess"
                )

        const_config_keys = ["params", "exec"]
        if not set(const_config_keys) == configd.keys():
            raise ConfigError(
                f"Config does not have correct top-level keys, expected {const_config_keys}, received {list(configd.keys())}"
            )

        try:
            excludes_opt = cast(bool, configd["exec"]["opts"]["mandate-excludes"])
        except KeyError:
            raise ConfigError("Config missing 'mandate-excludes' key from 'opts' section")
        if not isinstance(excludes_opt, bool):
            raise ConfigError(
                "Could not interpret 'mandate-excludes' key from 'opts' section as bool"
            )
        paramd = cast(dict[str, Any], configd["params"])
        execd = cast(dict[str, Any], configd["exec"])
        if not isinstance(execd, Mapping):
            raise ConfigError("Could not interpret 'exec' section - not a Mapping")
        if not isinstance(paramd, Mapping):
            raise ConfigError("Could not interpret 'params' value - not a Mapping")

        if not cls.check_namespacing_clashfree(proc_pool):  # pyright: ignore[reportArgumentType]
            raise ProcessError("Process pool with clashing namespaces")

        proc_nsd = {proc_type.ProcessNamespace: proc_type for proc_type in proc_pool}

        param_ns_set = set(paramd.keys())
        exec_ns_set = set(execd.keys())
        if unknown_proc_names := param_ns_set - set(proc_nsd.keys()):
            raise ConfigError(f"Params config requesting unknown processes: {unknown_proc_names}")
        if not param_ns_set.issubset(exec_ns_set):
            raise ConfigError(
                f"Param config keys are not a subset of exec config keys. Param only: {param_ns_set - exec_ns_set}"
            )

        ordered_proc_inst: list[ReadAwareProcessProtocol[Any, Any, FlagResult | None]] = []
        for proc_name in execd.keys():
            if proc_name == "opts":
                continue
            if not execd[proc_name]["enable"]:
                continue
            try:
                req_marks = execd[proc_name]["require-marks"]
            except KeyError:
                raise ConfigError(f"Config missing key 'require-marks' for process {proc_name}")
            try:
                exc_marks = execd[proc_name]["exclude-marks"]
            except KeyError:
                raise ConfigError(f"Config missing key 'exclude-marks' for process {proc_name}")

            proc_class = proc_nsd[proc_name]
            proc_fixed_param_class = proc_class.FixedParamClass
            proc_config_params = cast(dict[str, Any], paramd.get(proc_name, {}))  # TODO: isinstance

            if proc_fixed_param_class is not None:
                try:
                    proc_fixed_params = proc_fixed_param_class(**proc_config_params)
                except ValidationError as ex:
                    errl = ex.errors()
                    msg_nerr = f"{len(errl)} errors in params config for process {proc_name!r}:\n"
                    msg_locl = ""
                    for errd in errl:
                        submsg = f"{errd['loc'][0]!r} - {errd['msg']}\n"
                        msg_locl += submsg
                    msg = msg_nerr + msg_locl
                    raise ConfigError(msg)
            else:
                proc_fixed_params = None

            ordered_proc_inst.append(proc_class(proc_fixed_params, req_marks, exc_marks))

        _, (taggers, flaggers) = cls.validate_exec(
            ordered_proc_inst, excludes_opt, raise_on_fail=True
        )

        ra_ex = cls(taggers, flaggers, excludes_opt)
        return ra_ex

    def run(self, run_data: RunParams) -> tuple[FlagResult, ...]:
        # TODO: register a run params type at init and isinstance it
        for proc in self.taggers:
            proc(run_data)
            proc.reset()

        res: list[FlagResult] = []
        for proc in self.flaggers:
            res.append(proc(run_data))
            proc.reset()

        return tuple(res)
