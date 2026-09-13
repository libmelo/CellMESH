"""Shared numerical checks for observed and permuted scoring."""
from __future__ import annotations

import numpy as np
from scipy.stats import gmean


def _check_numeric_result(values, stage: str) -> None:
    """Reject invalid computed values before a mask can hide them."""
    array = np.asarray(values)
    if np.any(~np.isfinite(array)) or np.any(array < 0.0):
        raise ValueError(f"{stage} must be finite and non-negative")


def _reaction_activity(expression) -> np.ndarray:
    """Calculate the existing geometric mean, checking both sides of it."""
    _check_numeric_result(expression, "reaction pseudobulk expression")
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        activity = gmean(expression + 1.0, axis=1) - 1.0
    _check_numeric_result(activity, "reaction activity")
    return activity


def _positive_reference_scores(values, method: str, stage: str) -> tuple[np.ndarray, float]:
    """Normalize positive values without hiding overflowing intermediates."""
    x = np.asarray(values, dtype=float)
    _check_numeric_result(x, stage)
    score = np.zeros_like(x, dtype=float)
    positive = x > 0.0
    if not positive.any():
        # No positive reference is a defined absence, not a numerical failure.
        return score, np.nan

    selected = x[positive]
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        reference = float(np.mean(selected) if method == "mean" else np.median(selected))
    _check_numeric_result(reference, f"{stage} reference")
    with np.errstate(over="ignore", invalid="ignore"):
        denominator = selected + reference
    _check_numeric_result(denominator, f"{stage} normalization denominator")
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        score[positive] = selected / denominator
    _check_numeric_result(score, f"{stage} normalized scores")
    # 必须在除法、where/fillna/clip 或样本 nanmedian 掩盖异常前检查中间量。
    # 有限输入仍可能求和溢出；有限数 / Inf 会变成零，仅检查最终分数不够。
    # 观测与置换共用本函数，禁止为提速删掉这些检查或把异常转换为零。
    # 回归测试：tests/test_numeric_overflow.py。
    return score, reference
