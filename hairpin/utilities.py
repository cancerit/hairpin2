from typing import TextIO, Optional
import argparse
import re

def log_decision(
    msg: str,
    decision_lvl: int,
    log_lvl: int,
    log_file: Optional[TextIO]
) -> None:
    if log_file is not None and decision_lvl >= log_lvl:
        print(msg, file = log_file)


# Define a custom action to set a string and counter when the argument is present
class SetLogAndSeverity(argparse.Action):
        def __call__(self, parser, namespace, values, option_string=None):
            setattr(namespace, self.dest, values)

            # count ls to set sev
            if option_string is not None:
                pattern = re.compile('[^-l]')
                if pattern.search(option_string):
                    raise KeyError
                else:
                    sev = option_string.count('l')
            else:
                raise KeyError

            setattr(namespace, 'severity', sev)

