import numpy as np

# Bin edges
eta_bins = np.linspace(-2.5, 2.5, 11)
pt_bins = np.linspace(25.0, 100.0, 16)


# Input files
INPUT_FILES = {
    "MuonReco_nominal": "MuonReco_ScaleFactors_Nominal.txt",
    "MuonReco_systematic": "MuonReco_ScaleFactors_Systematic.txt",
    "MuonReco_statistic": "MuonReco_ScaleFactors_Statistic.txt",
    "MuonTrigger_nominal": "MuonTrigger_ScaleFactors_Nominal.txt",
    "MuonTrigger_systematic": "MuonTrigger_ScaleFactors_Systematic.txt",
    "MuonTrigger_statistic": "MuonTrigger_ScaleFactors_Statistic.txt",
    "MuonIso_nominal": "MuonIso_ScaleFactors_Nominal.txt",
    "MuonIso_systematic": "MuonIso_ScaleFactors_Systematic.txt",
    "MuonIso_statistic": "MuonIso_ScaleFactors_Statistic.txt",
}

# Internal storage (loaded on first use)
_MuonReco_nominal_sf = None
_MuonReco_systematic_sf = None
_MuonReco_statistic_sf = None
_MuonTrigger_nominal_sf = None
_MuonTrigger_systematic_sf = None
_MuonTrigger_statistic_sf = None
_MuonIso_nominal_sf = None
_MuonIso_systematic_sf = None
_MuonIso_statistic_sf = None


def _load_sf_maps():
    """Load the scale factor maps once."""
    global _MuonReco_nominal_sf, _MuonReco_systematic_sf, _MuonReco_statistic_sf, _MuonTrigger_nominal_sf, _MuonTrigger_systematic_sf, _MuonTrigger_statistic_sf, _MuonIso_nominal_sf, _MuonIso_systematic_sf, _MuonIso_statistic_sf

    if _MuonReco_nominal_sf is None:
        _MuonReco_nominal_sf = np.loadtxt(INPUT_FILES["MuonReco_nominal"], comments="#")
        _MuonReco_systematic_sf = np.loadtxt(INPUT_FILES["MuonReco_systematic"], comments="#")
        _MuonReco_statistic_sf = np.loadtxt(INPUT_FILES["MuonReco_statistic"], comments="#")
        _MuonTrigger_nominal_sf = np.loadtxt(INPUT_FILES["MuonTrigger_nominal"], comments="#")
        _MuonTrigger_systematic_sf = np.loadtxt(INPUT_FILES["MuonTrigger_systematic"], comments="#")
        _MuonTrigger_statistic_sf = np.loadtxt(INPUT_FILES["MuonTrigger_statistic"], comments="#")
        _MuonIso_nominal_sf = np.loadtxt(INPUT_FILES["MuonIso_nominal"], comments="#")
        _MuonIso_systematic_sf = np.loadtxt(INPUT_FILES["MuonIso_systematic"], comments="#")
        _MuonIso_statistic_sf = np.loadtxt(INPUT_FILES["MuonIso_statistic"], comments="#")


def getMuonRecoSF(pt, eta, variation=0):
    """
    Return the muon reconstruction scale factor.

    Parameters
    ----------
    pt : float
        Muon pT.
    eta : float
        Muon eta.
    variation : int, optional
        0 : nominal
        1 : nominal + systematic
        2 : nominal + statistical
        3 : nominal - systematic
        4 : nominal - statistical

    Returns
    -------
    float
        Requested scale factor.
    """

    _load_sf_maps()

    # Determine bin indices
    eta_bin = np.clip(np.digitize(eta, eta_bins) - 1, 0, len(eta_bins) - 2)
    pt_bin = np.clip(np.digitize(pt, pt_bins) - 1, 0, len(pt_bins) - 2)

    nominal = _MuonReco_nominal_sf[eta_bin, pt_bin]

    if variation == 0:
        return nominal
    elif variation == 1:
        return nominal + _MuonReco_systematic_sf[eta_bin, pt_bin]
    elif variation == 2:
        return nominal + _MuonReco_statistic_sf[eta_bin, pt_bin]
    elif variation == 3:
        return nominal - _MuonReco_systematic_sf[eta_bin, pt_bin]
    elif variation == 4:
        return nominal - _MuonReco_statistic_sf[eta_bin, pt_bin]
    else:
        raise ValueError(
            "variation must be one of {0, 1, 2, 3, 4}"
        )

def getMuonTriggerSF(pt, eta, variation=0):
    """
    Return the muon Trigger scale factor.
    """

    _load_sf_maps()

    # Determine bin indices
    eta_bin = np.clip(np.digitize(eta, eta_bins) - 1, 0, len(eta_bins) - 2)
    pt_bin = np.clip(np.digitize(pt, pt_bins) - 1, 0, len(pt_bins) - 2)

    nominal = _MuonTrigger_nominal_sf[eta_bin, pt_bin]

    if variation == 0:
        return nominal
    elif variation == 1:
        return nominal + _MuonTrigger_systematic_sf[eta_bin, pt_bin]
    elif variation == 2:
        return nominal + _MuonTrigger_statistic_sf[eta_bin, pt_bin]
    elif variation == 3:
        return nominal - _MuonTrigger_systematic_sf[eta_bin, pt_bin]
    elif variation == 4:
        return nominal - _MuonTrigger_statistic_sf[eta_bin, pt_bin]
    else:
        raise ValueError(
            "variation must be one of {0, 1, 2, 3, 4}"
        )


def getMuonIsolationSF(pt, eta, variation=0):
    """
    Return the muon Trigger scale factor.
    """

    _load_sf_maps()

    # Determine bin indices
    eta_bin = np.clip(np.digitize(eta, eta_bins) - 1, 0, len(eta_bins) - 2)
    pt_bin = np.clip(np.digitize(pt, pt_bins) - 1, 0, len(pt_bins) - 2)

    nominal = _MuonIso_nominal_sf[eta_bin, pt_bin]

    if variation == 0:
        return nominal
    elif variation == 1:
        return nominal + _MuonIso_systematic_sf[eta_bin, pt_bin]
    elif variation == 2:
        return nominal + _MuonIso_statistic_sf[eta_bin, pt_bin]
    elif variation == 3:
        return nominal - _MuonIso_systematic_sf[eta_bin, pt_bin]
    elif variation == 4:
        return nominal - _MuonIso_statistic_sf[eta_bin, pt_bin]
    else:
        raise ValueError(
            "variation must be one of {0, 1, 2, 3, 4}"
        )


if __name__ == "__main__":
    print("Muon Reco Scale Factor")
    sf = getMuonRecoSF(pt=32.3, eta=0.8)
    sf_up_sys = getMuonRecoSF(42.3, 0.8, variation=1)
    sf_down_stat = getMuonRecoSF(42.3, 0.8, variation=4)
    print(sf, sf_up_sys, sf_down_stat)

    print("Muon Trigger Scale Factor")
    sf = getMuonTriggerSF(pt=32.3, eta=0.8)
    sf_up_sys = getMuonTriggerSF(42.3, 0.8, variation=1)
    sf_down_stat = getMuonTriggerSF(42.3, 0.8, variation=4)
    print(sf, sf_up_sys, sf_down_stat)

    print("Muon Isolation Scale Factor")
    sf = getMuonIsolationSF(pt=32.3, eta=0.8)
    sf_up_sys = getMuonIsolationSF(42.3, 0.8, variation=1)
    sf_down_stat = getMuonIsolationSF(42.3, 0.8, variation=4)
    print(sf, sf_up_sys, sf_down_stat)