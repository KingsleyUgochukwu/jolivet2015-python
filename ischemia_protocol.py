
def ischemia_blood_flow(
    t,
    F0=0.012,
    onset=100.0,
    duration=1000.0,
    reduction=0.80,
):
    """
    Simulate transient cerebral ischemia and reperfusion.

    Parameters
    ----------
    t : float
        Simulation time in seconds.
    F0 : float
        Baseline cerebral blood flow.
    onset : float
        Time at which ischemia begins (seconds).
    duration : float
        Duration of ischemia (seconds).
    reduction : float
        Fractional reduction in cerebral blood flow.
        For example, 0.80 means an 80% reduction.

    Returns
    -------
    float
        Cerebral blood flow at time t.
    """

    if F0 <= 0:
        raise ValueError("Baseline blood flow must be positive.")

    if onset < 0 or duration <= 0:
        raise ValueError("Onset must be nonnegative and duration positive.")

    if not 0 <= reduction <= 1:
        raise ValueError("Reduction must be between 0 and 1.")

    if onset <= t < onset + duration:
        return F0 * (1.0 - reduction)

    return F0
