"""Shared numerical rules for observation, permutation and single-cell plots."""
from __future__ import annotations

from contextvars import ContextVar
from functools import wraps
import logging
import math
from threading import Lock

import numpy as np


_logger = logging.getLogger(__name__)
_underflow_reporter = ContextVar('cellmesh_underflow_reporter', default=None)


class _UnderflowReport:
    def __init__(self):
        self.stages = {}
        self.lock = Lock()

    def record(self, stage, count):
        with self.lock:
            self.stages[stage] = self.stages.get(stage, 0) + int(count)

    def snapshot(self):
        return {
            'underflow_detected': True,
            'policy': 'report_and_continue',
            'stages': [{'stage': stage, 'n_values': count}
                       for stage, count in sorted(self.stages.items())],
        }


def _report_underflow(mask, stage):
    """Report actual positive-to-zero loss; never threshold small positive values."""
    count = int(np.count_nonzero(mask))
    if not count:
        return
    report = _underflow_reporter.get()
    if report is not None:
        report.record(stage, count)
    else:
        # Logging is intentional: warnings configured as errors must not turn
        # the user's report-and-continue policy into an exception.
        _logger.warning('Numerical underflow in %s (%d values); continuing with '
                        'rounded zeros. These zeros do not establish absent expression.', stage, count)


def _with_numerical_diagnostics(function):
    """Aggregate one notice per public call, including nested/threaded scoring."""
    @wraps(function)
    def wrapped(*args, **kwargs):
        if _underflow_reporter.get() is not None:
            return function(*args, **kwargs)
        report = _UnderflowReport()
        token = _underflow_reporter.set(report)
        try:
            result = function(*args, **kwargs)
            if report.stages:
                diagnostic = report.snapshot()
                if hasattr(result, 'parameters'):
                    result.parameters['numerical_diagnostics'] = diagnostic
                    result.events.attrs['numerical_diagnostics'] = diagnostic
                elif isinstance(result, dict):
                    result['numerical_diagnostics'] = diagnostic
                elif hasattr(result, 'attrs'):
                    result.attrs['numerical_diagnostics'] = diagnostic
            return result
        finally:
            _underflow_reporter.reset(token)
            if report.stages:
                _logger.warning(
                    'Numerical underflow handled by continuing with rounded zeros. '
                    'Stages (value occurrences): %s. References, production states, '
                    'scores and p/FDR may be affected; zeros are not proof of absent '
                    'expression. See numerical_diagnostics in the returned result.',
                    ', '.join(f'{stage}={count}' for stage, count in sorted(report.stages.items())),
                )
    return wrapped


def _check_numeric_result(values, stage: str) -> None:
    """Reject non-finite/negative results; genuine underflow is reported separately."""
    array = np.asarray(values)
    if np.any(~np.isfinite(array)) or np.any(array < 0.0):
        raise ValueError(f"{stage} must be finite and non-negative")


def _positive_product(left, right, stage):
    left, right = np.asarray(left, dtype=float), np.asarray(right, dtype=float)
    with np.errstate(over='ignore', invalid='ignore', under='ignore'):
        result = left * right
    _check_numeric_result(result, stage)
    _report_underflow((left > 0) & (right > 0) & (result == 0), stage)
    return result


def _abundance_weights(fractions, exponent, stage='abundance weights'):
    fractions = np.asarray(fractions, dtype=float)
    with np.errstate(over='ignore', invalid='ignore', under='ignore'):
        weights = np.power(fractions, exponent)
    _check_numeric_result(weights, stage)
    _report_underflow((fractions > 0) & (weights == 0), stage)
    return weights


def _reaction_activity(expression, *, stage='reaction activity') -> np.ndarray:
    """Equal-weight gmean(1+x)-1, without losing tiny x when adding/subtracting 1."""
    expression = np.asarray(expression, dtype=float)
    _check_numeric_result(expression, 'reaction pseudobulk expression')
    # Do not restore gmean(x+1)-1 for speed: 1+1e-20 rounds to 1. Observation,
    # compiled permutations and cell plots must share this stable calculation.
    # A one-gene reaction is exactly x, also avoiding an unnecessary log roundtrip.
    with np.errstate(over='ignore', invalid='ignore', divide='ignore', under='ignore'):
        if expression.shape[1] == 1:
            activity = expression[:, 0].copy()
        else:
            activity = np.expm1(np.mean(np.log1p(expression), axis=1))
    _check_numeric_result(activity, 'reaction activity')
    _report_underflow(np.any(expression > 0, axis=1) & (activity == 0), stage)
    return activity


def _sender_base(production, consumption, stage='sender normalization'):
    production, consumption = np.asarray(production, dtype=float), np.asarray(consumption, dtype=float)
    with np.errstate(over='ignore', invalid='ignore', under='ignore'):
        denominator = production + consumption
    _check_numeric_result(denominator, f'{stage} denominator')
    with np.errstate(over='ignore', invalid='ignore', divide='ignore', under='ignore'):
        ratio = np.divide(production, denominator, out=np.zeros_like(production), where=denominator > 0)
    _check_numeric_result(ratio, f'{stage} ratio')
    _report_underflow((production > 0) & (denominator > 0) & (ratio == 0), f'{stage} ratio')
    # Squaring P first loses representable results, e.g. P=2e-200, C=0.
    return _positive_product(production, ratio, f'{stage} base')


def _event_score(sender, receiver, stage='event scores'):
    # Event-table construction calls this for scalar pairs. Keep scalar checks
    # inexpensive without restoring sqrt(sender*receiver), which loses small
    # representable scores. Both paths still use the same double-precision roots.
    if isinstance(sender, (int, float, np.number)) and isinstance(receiver, (int, float, np.number)):
        if not (math.isfinite(sender) and math.isfinite(receiver)) or sender < 0 or receiver < 0:
            raise ValueError(f'{stage} must be finite and non-negative')
        result = math.sqrt(sender) * math.sqrt(receiver)
        if not math.isfinite(result):
            raise ValueError(f'{stage} must be finite and non-negative')
        if sender > 0 and receiver > 0 and result == 0:
            _report_underflow(True, stage)
        return result
    _check_numeric_result(sender, stage)
    _check_numeric_result(receiver, stage)
    with np.errstate(under='ignore'):
        # The product can underflow even though its square root is representable.
        return _positive_product(np.sqrt(sender), np.sqrt(receiver), stage)


def _positive_reference_scores(values, method: str, stage: str) -> tuple[np.ndarray, float]:
    """Normalize positive values without hiding overflowing intermediates."""
    x = np.asarray(values, dtype=float)
    _check_numeric_result(x, stage)
    score = np.zeros_like(x, dtype=float)
    positive = x > 0.0
    if not positive.any():
        return score, np.nan
    selected = x[positive]
    with np.errstate(over='ignore', invalid='ignore', divide='ignore', under='ignore'):
        reference = float(np.mean(selected) if method == 'mean' else np.median(selected))
    _check_numeric_result(reference, f'{stage} reference')
    _report_underflow(reference == 0, f'{stage} reference')
    with np.errstate(over='ignore', invalid='ignore', under='ignore'):
        denominator = selected + reference
    _check_numeric_result(denominator, f'{stage} normalization denominator')
    with np.errstate(over='ignore', invalid='ignore', divide='ignore', under='ignore'):
        score[positive] = selected / denominator
    _check_numeric_result(score, f'{stage} normalized scores')
    _report_underflow(score[positive] == 0, f'{stage} normalized scores')
    # Keep strict-positive references and the overflow checks in both paths.
    # Do not replace rounding zeros by an arbitrary epsilon or skip permutations.
    # Regression coverage: test_numeric_stability.py and test_numeric_overflow.py.
    return score, reference
