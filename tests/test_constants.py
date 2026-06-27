from openfep.constants import (
    FORCEFIELD,
    StrucTag,
    KEEP_STRUC_TAGS,
    COMPLEX_BUFFER_WIDTH,
    NET_CHARGE_COMPLEX_BUFFER_WIDTH,
    SOLVENT_BUFFER_WIDTH,
    VACUUM_BUFFER_WIDTH,
    HMR_TIMESTEPS,
    HMR_MIGRATION_INTERVAL,
    MIN_CHARGED_SALT_CONC,
    GCMC_SOLVENT_VDW_SCALE_FACTOR,
    ENEGRP_NAME_TEMPLATE,
    ENEGRP_FIRST,
    ENEGRP_INTERVAL,
    ABFE_COMPLEX_LAMBDAS,
    ABFE_SOLVENT_LAMBDAS,
    ABFE_RESTRAINED_COMPLEX_LAMBDAS,
    ABFE_RESTRAINED_SOLVENT_LAMBDAS,
    ABFE_CHARGED_COMPLEX_LAMBDAS,
    ABFE_CHARGED_SOLVENT_LAMBDAS,
    ABFE_CHARGED_RESTRAINED_COMPLEX_LAMBDAS,
    ABFE_CHARGED_RESTRAINED_SOLVENT_LAMBDAS,
    RBFE_DEFAULT_LAMBDAS,
    BACKBONE_ASL,
    RECEPTOR_ASL,
)


def test_forcefield_is_opls4():
    assert FORCEFIELD == "OPLS4"


def test_struc_tag_values():
    assert StrucTag.RECEPTOR == "receptor"
    assert StrucTag.LIGAND == "ligand"
    assert StrucTag.COMPLEX == "complex"
    assert StrucTag.SOLVENT == "solvent"
    assert StrucTag.MEMBRANE == "membrane"


def test_keep_struc_tags_rbfe_complex_includes_membrane():
    tags = KEEP_STRUC_TAGS["rbfe"]["complex"]
    assert StrucTag.MEMBRANE in tags
    assert StrucTag.RECEPTOR in tags
    assert StrucTag.LIGAND in tags


def test_keep_struc_tags_abfe_md_includes_membrane():
    tags = KEEP_STRUC_TAGS["abfe"]["md"]
    assert StrucTag.MEMBRANE in tags
    assert StrucTag.RECEPTOR in tags
    assert StrucTag.LIGAND in tags


def test_abfe_complex_leg_keeps_complex_not_ligand():
    tags = KEEP_STRUC_TAGS["abfe"]["complex"]
    assert StrucTag.COMPLEX in tags
    assert StrucTag.LIGAND not in tags


def test_buffer_widths():
    assert COMPLEX_BUFFER_WIDTH == 5.0
    assert NET_CHARGE_COMPLEX_BUFFER_WIDTH == 8.0
    assert SOLVENT_BUFFER_WIDTH == 10.0
    assert VACUUM_BUFFER_WIDTH == 100.0


def test_hmr():
    assert HMR_TIMESTEPS == [0.004, 0.004, 0.008]
    assert HMR_MIGRATION_INTERVAL == 0.024


def test_abfe_lambda_counts_default():
    assert ABFE_COMPLEX_LAMBDAS == 68
    assert ABFE_SOLVENT_LAMBDAS == 60
    assert ABFE_RESTRAINED_COMPLEX_LAMBDAS == 80
    assert ABFE_RESTRAINED_SOLVENT_LAMBDAS == 68


def test_abfe_lambda_counts_charged():
    assert ABFE_CHARGED_COMPLEX_LAMBDAS == 108
    assert ABFE_CHARGED_SOLVENT_LAMBDAS == 60
    assert ABFE_CHARGED_RESTRAINED_COMPLEX_LAMBDAS == 128
    assert ABFE_CHARGED_RESTRAINED_SOLVENT_LAMBDAS == 68


def test_backbone_asl_includes_nucleic():
    assert "nucleic_acids" in BACKBONE_ASL
    assert "protein" in BACKBONE_ASL
    assert "not a.ele H" in BACKBONE_ASL


def test_enegrp_template():
    assert "$JOBNAME$" in ENEGRP_NAME_TEMPLATE
    assert "_enegrp.dat" in ENEGRP_NAME_TEMPLATE
