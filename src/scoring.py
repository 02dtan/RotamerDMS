"""
Multi-objective scoring and Pareto front analysis for rotamer optimization.

Provides:
- Z-score normalization for combining metrics with different scales
- Joint scoring function combining volume, deltaG, and WCA
- Pareto front computation for multi-objective tradeoff analysis
"""

import math
from typing import List, Dict, Tuple, Optional


def compute_stats(values: List[float]) -> Tuple[float, float]:
    """
    Compute mean and standard deviation for z-score normalization.
    
    Args:
        values: List of numeric values (None values are skipped)
        
    Returns:
        Tuple of (mean, std) for use in z-score computation
    """
    valid_values = [v for v in values if v is not None]
    
    if len(valid_values) < 2:
        # Not enough data for meaningful statistics
        return 0.0, 1.0
    
    mean = sum(valid_values) / len(valid_values)
    variance = sum((v - mean) ** 2 for v in valid_values) / len(valid_values)
    std = math.sqrt(variance) if variance > 0 else 1.0
    
    return mean, std


def compute_joint_score(
    volume_change: float,
    delta_g: Optional[float],
    wca_energy: Optional[float],
    w_vol: float = 1.0,
    w_dg: float = 1.0,
    w_wca: float = 1.0,
    vol_mean: float = 0.0,
    vol_std: float = 1.0,
    dg_mean: float = 0.0,
    dg_std: float = 1.0,
    wca_mean: float = 0.0,
    wca_std: float = 1.0,
    wca_threshold: Optional[float] = None
) -> Optional[float]:
    """
    Compute joint optimization score for a rotamer state.
    
    Score = w_vol * z(ΔV) - w_dg * z(ΔG) - w_wca * z(WCA)
    
    Higher score is better (more volume, less deltaG, less WCA).
    
    Args:
        volume_change: Volume change in Å³
        delta_g: Conformational free energy change (arbitrary units)
        wca_energy: WCA repulsive energy (REU)
        w_vol: Weight for volume term (default 1.0)
        w_dg: Weight for deltaG term (default 1.0)
        w_wca: Weight for WCA term (default 1.0)
        vol_mean, vol_std: Normalization parameters for volume
        dg_mean, dg_std: Normalization parameters for deltaG
        wca_mean, wca_std: Normalization parameters for WCA
        wca_threshold: If set, return None for rotamers with WCA > threshold (hard clash filter)
        
    Returns:
        Joint score (higher is better), or None if hard clash threshold exceeded
    """
    # Hard clash filter
    if wca_threshold is not None and wca_energy is not None:
        if wca_energy > wca_threshold:
            return None
    
    # Normalize to z-scores
    z_vol = (volume_change - vol_mean) / vol_std if vol_std > 0 else 0.0
    
    z_dg = 0.0
    if delta_g is not None and dg_std > 0:
        z_dg = (delta_g - dg_mean) / dg_std
    
    z_wca = 0.0
    if wca_energy is not None and wca_std > 0:
        z_wca = (wca_energy - wca_mean) / wca_std
    
    # Combined score: maximize volume, minimize deltaG and WCA
    score = w_vol * z_vol - w_dg * z_dg - w_wca * z_wca
    
    return score


def compute_normalization_params(rotamer_results: List[Dict]) -> Dict:
    """
    Compute normalization parameters (mean, std) for each metric.
    
    Args:
        rotamer_results: List of rotamer result dictionaries
        
    Returns:
        Dictionary with normalization parameters for each metric
    """
    volumes = [r['volume_change'] for r in rotamer_results]
    delta_gs = [r['delta_g'] for r in rotamer_results]
    wcas = [r['wca_energy'] for r in rotamer_results]
    
    vol_mean, vol_std = compute_stats(volumes)
    dg_mean, dg_std = compute_stats(delta_gs)
    wca_mean, wca_std = compute_stats(wcas)
    
    return {
        'vol_mean': vol_mean, 'vol_std': vol_std,
        'dg_mean': dg_mean, 'dg_std': dg_std,
        'wca_mean': wca_mean, 'wca_std': wca_std
    }


def rank_rotamers_by_joint_score(
    rotamer_results: List[Dict],
    w_vol: float = 1.0,
    w_dg: float = 1.0,
    w_wca: float = 1.0,
    wca_threshold: Optional[float] = None
) -> List[Dict]:
    """
    Rank rotamers by joint optimization score.
    
    Args:
        rotamer_results: List of rotamer result dictionaries
        w_vol: Weight for volume (default 1.0)
        w_dg: Weight for deltaG (default 1.0)
        w_wca: Weight for WCA (default 1.0)
        wca_threshold: Hard clash filter threshold (REU)
        
    Returns:
        List of rotamer results sorted by joint score (descending),
        with 'joint_score' field added to each
    """
    if not rotamer_results:
        return []
    
    # Compute normalization parameters
    norm_params = compute_normalization_params(rotamer_results)
    
    # Add joint scores
    scored_results = []
    for r in rotamer_results:
        score = compute_joint_score(
            r['volume_change'],
            r['delta_g'],
            r['wca_energy'],
            w_vol=w_vol, w_dg=w_dg, w_wca=w_wca,
            wca_threshold=wca_threshold,
            **norm_params
        )
        
        result_copy = r.copy()
        result_copy['joint_score'] = score
        scored_results.append(result_copy)
    
    # Filter out None scores (hard clash rejects) and sort
    valid_results = [r for r in scored_results if r['joint_score'] is not None]
    valid_results.sort(key=lambda x: x['joint_score'], reverse=True)
    
    # Append rejected results at the end (marked with None score)
    rejected = [r for r in scored_results if r['joint_score'] is None]
    
    return valid_results + rejected


def compute_pareto_front_2d(
    rotamer_results: List[Dict],
    maximize_key: str = 'volume_change',
    minimize_key: str = 'wca_energy'
) -> List[Dict]:
    """
    Compute 2D Pareto front for volume vs. WCA tradeoff.
    
    A rotamer is Pareto-optimal if no other rotamer is better in BOTH objectives.
    
    Args:
        rotamer_results: List of rotamer result dictionaries
        maximize_key: Key to maximize (default 'volume_change')
        minimize_key: Key to minimize (default 'wca_energy')
        
    Returns:
        List of Pareto-optimal rotamer results, sorted by maximize_key descending
    """
    # Filter to results with valid values for both keys
    valid_results = [
        r for r in rotamer_results 
        if r.get(maximize_key) is not None and r.get(minimize_key) is not None
    ]
    
    if not valid_results:
        return []
    
    # Sort by maximize_key descending, then by minimize_key ascending (for ties)
    # This ensures that among equal maximize_key values, the best minimize_key comes first
    sorted_results = sorted(valid_results, key=lambda x: (-x[maximize_key], x[minimize_key]))
    
    # Sweep to find Pareto front
    pareto_front = []
    min_seen = float('inf')  # Best (minimum) value of minimize_key seen so far
    
    for r in sorted_results:
        # A point is on the Pareto front if its minimize_key value is better
        # than any point with higher maximize_key value
        if r[minimize_key] < min_seen:
            pareto_front.append(r)
            min_seen = r[minimize_key]
    
    return pareto_front


def compute_pareto_front_3d(
    rotamer_results: List[Dict],
    maximize_keys: List[str] = None,
    minimize_keys: List[str] = None
) -> List[Dict]:
    """
    Compute 3D Pareto front for multi-objective optimization.
    
    Args:
        rotamer_results: List of rotamer result dictionaries
        maximize_keys: Keys to maximize (default ['volume_change'])
        minimize_keys: Keys to minimize (default ['delta_g', 'wca_energy'])
        
    Returns:
        List of Pareto-optimal rotamer results
    """
    if maximize_keys is None:
        maximize_keys = ['volume_change']
    if minimize_keys is None:
        minimize_keys = ['delta_g', 'wca_energy']
    
    all_keys = maximize_keys + minimize_keys
    
    # Filter to results with valid values for all keys
    valid_results = [
        r for r in rotamer_results 
        if all(r.get(k) is not None for k in all_keys)
    ]
    
    if not valid_results:
        return []
    
    def dominates(a: Dict, b: Dict) -> bool:
        """Check if solution a dominates solution b (a is at least as good in all, strictly better in one)."""
        a_at_least_as_good_in_all = True
        a_strictly_better_in_one = False
        
        for key in maximize_keys:
            if a[key] < b[key]:  # a is worse
                a_at_least_as_good_in_all = False
            if a[key] > b[key]:  # a is better
                a_strictly_better_in_one = True
        
        for key in minimize_keys:
            if a[key] > b[key]:  # a is worse (higher = worse for minimize)
                a_at_least_as_good_in_all = False
            if a[key] < b[key]:  # a is better
                a_strictly_better_in_one = True
        
        return a_at_least_as_good_in_all and a_strictly_better_in_one
    
    # Find non-dominated solutions
    pareto_front = []
    for candidate in valid_results:
        is_dominated = False
        for other in valid_results:
            if candidate is not other and dominates(other, candidate):
                is_dominated = True
                break
        if not is_dominated:
            pareto_front.append(candidate)
    
    # Sort by first maximize_key descending
    pareto_front.sort(key=lambda x: x[maximize_keys[0]], reverse=True)
    
    return pareto_front


def analyze_pareto_tradeoff(
    rotamer_results: List[Dict],
    verbose: bool = True
) -> Dict:
    """
    Comprehensive Pareto analysis for rotamer selection.
    
    Args:
        rotamer_results: List of rotamer result dictionaries
        verbose: Print analysis summary
        
    Returns:
        Dictionary containing:
            - pareto_2d_vol_wca: 2D Pareto front (volume vs WCA)
            - pareto_3d: 3D Pareto front (volume vs deltaG vs WCA)
            - best_volume: Rotamer with best volume
            - best_wca: Rotamer with lowest WCA
            - best_balanced: Rotamer on Pareto front closest to ideal
    """
    # 2D Pareto: volume vs WCA
    pareto_2d = compute_pareto_front_2d(rotamer_results, 'volume_change', 'wca_energy')
    
    # 3D Pareto: volume vs deltaG vs WCA
    pareto_3d = compute_pareto_front_3d(
        rotamer_results,
        maximize_keys=['volume_change'],
        minimize_keys=['delta_g', 'wca_energy']
    )
    
    # Find extremes
    valid_vol = [r for r in rotamer_results if r.get('volume_change') is not None]
    valid_wca = [r for r in rotamer_results if r.get('wca_energy') is not None]
    
    best_volume = max(valid_vol, key=lambda x: x['volume_change']) if valid_vol else None
    best_wca = min(valid_wca, key=lambda x: x['wca_energy']) if valid_wca else None
    
    # Best balanced: Pareto point with best joint score (equal weights)
    if pareto_2d:
        scored_pareto = rank_rotamers_by_joint_score(pareto_2d)
        best_balanced = scored_pareto[0] if scored_pareto else None
    else:
        best_balanced = None
    
    if verbose and pareto_2d:
        print(f"\n    Pareto Front Analysis (Volume vs WCA):")
        print(f"    {len(pareto_2d)} Pareto-optimal rotamers out of {len(rotamer_results)}")
        print(f"    Pareto front (volume ↓, WCA →):")
        for i, r in enumerate(pareto_2d[:5]):  # Show top 5
            print(f"      {i+1}. idx={r['index']}: ΔV={r['volume_change']:+.2f} Å³, "
                  f"WCA={r['wca_energy']:.2f} REU")
        if len(pareto_2d) > 5:
            print(f"      ... and {len(pareto_2d) - 5} more")
    
    return {
        'pareto_2d_vol_wca': pareto_2d,
        'pareto_3d': pareto_3d,
        'best_volume': best_volume,
        'best_wca': best_wca,
        'best_balanced': best_balanced,
        'num_pareto_2d': len(pareto_2d),
        'num_pareto_3d': len(pareto_3d)
    }
