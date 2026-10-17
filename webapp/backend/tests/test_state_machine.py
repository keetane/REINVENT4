import pytest

from app.models.job import ALLOWED_TRANSITIONS, JobStatus, TransitionError, require_transition


@pytest.mark.parametrize(
    ("current", "target"),
    [
        (current, target)
        for current, targets in ALLOWED_TRANSITIONS.items()
        for target in targets
    ],
)
def test_all_declared_transitions_are_allowed(current, target):
    require_transition(current, target)


@pytest.mark.parametrize(
    ("current", "target"),
    [
        (current, target)
        for current in JobStatus
        for target in JobStatus
        if target not in ALLOWED_TRANSITIONS[current]
    ],
)
def test_all_other_transitions_are_rejected(current, target):
    with pytest.raises(TransitionError):
        require_transition(current, target)

