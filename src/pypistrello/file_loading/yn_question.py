#
# Copyright 2026 Universidad Complutense de Madrid
#
# This file is part of pypistrello.
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE
#


def question_yes_no(question, default=False):
    """
    Ask a yes/no question via input() and return True or False.

    Parameters
    ----------
    question : str
        Question to display.
    default : bool
        Default answer if user presses Enter.

    Returns
    -------
    bool
        True for yes, False for no.
    """
    prompt = " (y/[n]): " if not default else " ([y]/n): " # Default: press Enter-> False and y-> True. If not default, the opposite

    while True:
        answer = input(question + prompt).strip().lower()

        if answer == "":
            return default
        elif answer in ("y", "yes"):
            return True
        elif answer in ("n", "no"):
            return False
        else:
            print("Please answer with 'y' or 'n'.")

