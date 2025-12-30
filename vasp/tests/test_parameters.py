"""Tests for parameter preset functions."""

import pytest

from vasp.parameters import (
    COMMON_U_VALUES,
    VDW_METHODS,
    HubbardU,
    get_bse_params,
    get_gw_params,
    get_hybrid_params,
    get_ldau_params,
    get_md_params,
    get_mlff_params,
    get_optical_params,
    get_phonon_params,
    get_soc_params,
    get_vdw_params,
)


class TestVdWParameters:
    """Test van der Waals correction parameters."""

    def test_vdw_methods_dict(self):
        """Test VDW_METHODS dictionary contains expected methods."""
        assert "d2" in VDW_METHODS
        assert "d3" in VDW_METHODS
        assert "d3bj" in VDW_METHODS
        assert "ts" in VDW_METHODS
        assert "vdw-df" in VDW_METHODS
        assert "vdw-df2" in VDW_METHODS

    def test_get_vdw_params_d3bj(self):
        """Test DFT-D3-BJ parameters."""
        params = get_vdw_params("d3bj")
        assert params["ivdw"] == 12
        assert "luse_vdw" not in params

    def test_get_vdw_params_d2(self):
        """Test DFT-D2 parameters."""
        params = get_vdw_params("d2")
        assert params["ivdw"] == 10  # Grimme D2

    def test_get_vdw_params_ts(self):
        """Test Tkatchenko-Scheffler parameters."""
        params = get_vdw_params("ts")
        assert params["ivdw"] == 2

    def test_get_vdw_params_vdw_df(self):
        """Test VdW-DF parameters."""
        params = get_vdw_params("vdw-df")
        assert params["luse_vdw"] is True
        assert params["aggac"] == 0.0
        assert "ivdw" not in params

    def test_get_vdw_params_vdw_df2(self):
        """Test VdW-DF2 parameters."""
        params = get_vdw_params("vdw-df2")
        assert params["luse_vdw"] is True
        assert params["gga"] == "ML"

    def test_get_vdw_params_invalid(self):
        """Test invalid VdW method raises error."""
        with pytest.raises(ValueError, match="Unknown vdW method"):
            get_vdw_params("invalid_method")

    def test_get_vdw_params_case_insensitive(self):
        """Test method names are case-insensitive."""
        params_lower = get_vdw_params("d3bj")
        params_upper = get_vdw_params("D3BJ")
        assert params_lower == params_upper


class TestLDAUParameters:
    """Test DFT+U parameters."""

    def test_common_u_values(self):
        """Test COMMON_U_VALUES contains transition metals."""
        assert "Fe" in COMMON_U_VALUES
        assert "Co" in COMMON_U_VALUES
        assert "Ni" in COMMON_U_VALUES
        assert "Ti" in COMMON_U_VALUES
        assert "V" in COMMON_U_VALUES
        assert "Mn" in COMMON_U_VALUES

    def test_common_u_values_lanthanides(self):
        """Test COMMON_U_VALUES contains lanthanides."""
        assert "Ce" in COMMON_U_VALUES
        assert "Pr" in COMMON_U_VALUES
        assert "Nd" in COMMON_U_VALUES

    def test_hubbardu_dataclass(self):
        """Test HubbardU dataclass."""
        hu = HubbardU(u=4.0, j=0.0, l_angular=2)
        assert hu.u == 4.0
        assert hu.j == 0.0
        assert hu.l_angular == 2

    def test_hubbardu_defaults(self):
        """Test HubbardU default values."""
        hu = HubbardU(u=3.5)
        assert hu.u == 3.5
        assert hu.j == 0.0
        assert hu.l_angular == 2  # Default to d-orbitals

    def test_get_ldau_params_single_element(self):
        """Test LDAU parameters for single element."""
        params = get_ldau_params(["Fe", "O"], {"Fe": HubbardU(u=4.0)})
        assert params["ldau"] is True
        assert params["ldautype"] == 2
        # O is in COMMON_U_VALUES with l=1 (p-orbital), u=0
        assert params["ldaul"] == [2, 1]  # Fe has d, O has p from COMMON_U_VALUES
        assert params["ldauu"] == [4.0, 0.0]
        assert params["ldauj"] == [0.0, 0.0]

    def test_get_ldau_params_multiple_elements(self):
        """Test LDAU parameters for multiple elements."""
        params = get_ldau_params(
            ["Fe", "Ni", "O"],
            {"Fe": HubbardU(u=4.0), "Ni": HubbardU(u=6.0, j=1.0)},
        )
        assert params["ldau"] is True
        assert params["ldaul"] == [2, 2, 1]  # O has l=1 from COMMON_U_VALUES
        assert params["ldauu"] == [4.0, 6.0, 0.0]
        assert params["ldauj"] == [0.0, 1.0, 0.0]

    def test_get_ldau_params_f_orbitals(self):
        """Test LDAU parameters for f-orbital elements."""
        params = get_ldau_params(["Ce", "O"], {"Ce": HubbardU(u=5.0, l_angular=3)})
        assert params["ldaul"] == [3, 1]  # Ce has f, O has l=1 from COMMON_U_VALUES

    def test_get_ldau_params_custom_type(self):
        """Test LDAU with custom LDAUTYPE."""
        params = get_ldau_params(
            ["Fe", "O"],
            {"Fe": HubbardU(u=4.0)},
            ldautype=1,  # Liechtenstein
        )
        assert params["ldautype"] == 1


class TestHybridParameters:
    """Test hybrid functional parameters."""

    def test_get_hybrid_params_hse06(self):
        """Test HSE06 parameters."""
        params = get_hybrid_params("hse06")
        assert params["lhfcalc"] is True
        assert params["hfscreen"] == 0.2
        assert params["algo"] == "Damped"  # Recommended for hybrids
        assert "precfock" in params

    def test_get_hybrid_params_pbe0(self):
        """Test PBE0 parameters."""
        params = get_hybrid_params("pbe0")
        assert params["lhfcalc"] is True
        assert params["aexx"] == 0.25
        assert params["algo"] == "Damped"
        # PBE0 has no screening (hfscreen not set or 0)
        assert params.get("hfscreen", 0) == 0

    def test_get_hybrid_params_b3lyp(self):
        """Test B3LYP parameters."""
        params = get_hybrid_params("b3lyp")
        assert params["lhfcalc"] is True
        assert params["aexx"] == 0.2
        assert "aggax" in params
        assert "aggac" in params

    def test_get_hybrid_params_kpts_reduction(self):
        """Test hybrid with k-point reduction."""
        params = get_hybrid_params("hse06", nkpts_factor=2)
        assert params["nkred"] == 2

    def test_get_hybrid_params_invalid(self):
        """Test invalid hybrid functional raises error."""
        with pytest.raises(ValueError, match="Unknown hybrid functional"):
            get_hybrid_params("invalid_hybrid")


class TestSOCParameters:
    """Test spin-orbit coupling parameters."""

    def test_get_soc_params_basic(self):
        """Test basic SOC parameters."""
        params = get_soc_params()
        assert params["lsorbit"] is True
        assert params["lnoncollinear"] is True
        assert params["gga_compat"] is True  # Default in implementation
        assert params["lmaxmix"] == 4

    def test_get_soc_params_with_saxis(self):
        """Test SOC with spin quantization axis."""
        params = get_soc_params(saxis=(0, 0, 1))
        assert params["saxis"] == [0, 0, 1]

    def test_get_soc_params_lorbmom(self):
        """Test SOC with orbital moment calculation."""
        params = get_soc_params(lorbmom=True)
        assert params["lorbmom"] is True


class TestMLFFParameters:
    """Test machine learning force field parameters."""

    def test_get_mlff_params_train(self):
        """Test MLFF training parameters."""
        params = get_mlff_params(mode="train")
        assert params["ml_lmlff"] is True
        assert params["ml_mode"] == "run"  # 'train' maps to 'run' in VASP
        assert "ml_rcut1" in params
        assert "ml_rcut2" in params
        # Training mode should have weights
        assert "ml_wtifor" in params

    def test_get_mlff_params_select(self):
        """Test MLFF selection parameters."""
        params = get_mlff_params(mode="select")
        assert params["ml_mode"] == "select"

    def test_get_mlff_params_run(self):
        """Test MLFF run parameters."""
        params = get_mlff_params(mode="run")
        assert params["ml_mode"] == "run"
        assert params["ml_lmlff"] is True

    def test_get_mlff_params_refit(self):
        """Test MLFF refit parameters."""
        params = get_mlff_params(mode="refit")
        assert params["ml_mode"] == "refit"

    def test_get_mlff_params_custom_cutoffs(self):
        """Test MLFF with custom cutoffs."""
        params = get_mlff_params(rcut1=10.0, rcut2=8.0)
        assert params["ml_rcut1"] == 10.0
        assert params["ml_rcut2"] == 8.0

    def test_get_mlff_params_invalid_mode(self):
        """Test invalid MLFF mode raises error."""
        with pytest.raises(ValueError, match="Unknown MLFF mode"):
            get_mlff_params(mode="invalid")


class TestMDParameters:
    """Test molecular dynamics parameters."""

    def test_get_md_params_nve(self):
        """Test NVE ensemble parameters."""
        params = get_md_params("nve", temperature=300.0, timestep=1.0, nsteps=1000)
        assert params["ibrion"] == 0
        assert params["mdalgo"] == 0
        assert params["tebeg"] == 300.0
        assert params["potim"] == 1.0
        assert params["nsw"] == 1000
        assert params["smass"] == -3  # NVE microcanonical

    def test_get_md_params_nvt_nose(self):
        """Test NVT Nose-Hoover parameters."""
        params = get_md_params("nvt-nose", temperature=500.0, timestep=2.0)
        assert params["ibrion"] == 0
        assert params["mdalgo"] == 2
        assert params["tebeg"] == 500.0
        assert params["smass"] == 0  # Nose-Hoover thermostat

    def test_get_md_params_nvt_langevin(self):
        """Test NVT Langevin parameters."""
        params = get_md_params("nvt-langevin", temperature=300.0)
        assert params["mdalgo"] == 3
        # Langevin doesn't require langevin_gamma in base params

    def test_get_md_params_npt(self):
        """Test NPT ensemble parameters."""
        params = get_md_params("npt-parrinello", temperature=300.0)
        assert params["mdalgo"] == 3
        assert params["isif"] == 3
        # Parrinello-Rahman NPT

    def test_get_md_params_invalid(self):
        """Test invalid ensemble raises error."""
        with pytest.raises(ValueError, match="Unknown MD ensemble"):
            get_md_params("invalid_ensemble")


class TestPhononParameters:
    """Test phonon calculation parameters."""

    def test_get_phonon_params_dfpt(self):
        """Test DFPT phonon parameters."""
        params = get_phonon_params("dfpt")
        assert params["ibrion"] == 8
        assert params["lepsilon"] is True  # Born charges
        assert params["addgrid"] is True

    def test_get_phonon_params_phonopy(self):
        """Test Phonopy phonon parameters."""
        params = get_phonon_params("phonopy")
        assert params["ibrion"] == -1  # Static for force calculation
        assert params["nsw"] == 0
        assert params["isym"] == 0

    def test_get_phonon_params_finite_diff(self):
        """Test finite difference phonon parameters."""
        params = get_phonon_params("finite-diff")
        assert params["ibrion"] == 6  # Central differences
        assert params["potim"] < 0.1  # Small displacement
        assert params["nfree"] == 2

    def test_get_phonon_params_invalid(self):
        """Test invalid phonon method raises error."""
        with pytest.raises(ValueError, match="Unknown phonon method"):
            get_phonon_params("invalid_method")


class TestOpticalParameters:
    """Test optical property parameters."""

    def test_get_optical_params(self):
        """Test optical calculation parameters."""
        params = get_optical_params()
        assert params["loptics"] is True
        assert "nedos" in params
        assert params["nedos"] > 1000  # Need many points for optics


class TestGWParameters:
    """Test GW approximation parameters."""

    def test_get_gw_params_g0w0(self):
        """Test G0W0 parameters (same as gw0)."""
        params = get_gw_params("gw0")
        assert params["algo"] == "GW0"
        assert "nomega" in params
        assert params["loptics"] is True

    def test_get_gw_params_gw(self):
        """Test full GW parameters."""
        params = get_gw_params("gw")
        assert params["algo"] == "GW"
        assert "nomega" in params

    def test_get_gw_params_scgw(self):
        """Test self-consistent GW parameters."""
        params = get_gw_params("scgw")
        assert params["algo"] == "scGW"

    def test_get_gw_params_invalid(self):
        """Test invalid GW method raises error."""
        with pytest.raises(ValueError, match="Unknown GW algorithm"):
            get_gw_params("invalid_gw")


class TestBSEParameters:
    """Test BSE (Bethe-Salpeter Equation) parameters."""

    def test_get_bse_params(self):
        """Test BSE parameters."""
        params = get_bse_params()
        assert params["algo"] == "BSE"
        assert params["antires"] == 0

    def test_get_bse_params_with_bands(self):
        """Test BSE with explicit band counts."""
        params = get_bse_params(nbandso=10, nbandsv=20)
        assert params["nbandso"] == 10
        assert params["nbandsv"] == 20
