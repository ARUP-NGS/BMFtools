
"""
Contains custom exceptions and other tools for handling errors,
particularly with cython's unique error-handling tools.
"""


class PermissionException(Exception):
    """
    Thrown when permissions are likely the preventing culprit
     for a given operation.
    """
    pass


class IllegalArgumentError(ValueError):
    pass


class ThisIsMadness(Exception):
    """
    Thrown when something just doesn't seem right and existing exceptions
    don't have enough panache.
    """
    pass