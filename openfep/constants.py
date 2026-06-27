# Forcefield
FORCEFIELD = "OPLS4"  # OPLSVersion.F17 string in MSJ. Not OPLS3e (F16).

# Buffer widths (Å) — from desmond/msj_constants.py
COMPLEX_BUFFER_WIDTH = 5.0
NET_CHARGE_COMPLEX_BUFFER_WIDTH = 8.0
SOLVENT_BUFFER_WIDTH = 10.0
VACUUM_BUFFER_WIDTH = 100.0

# HMR — from desmond/msj_constants.py _HMR_TIMESTEPS / _HMR_MIGRATION_INTERVAL
HMR_TIMESTEPS = [0.004, 0.004, 0.008]   # ps
HMR_MIGRATION_INTERVAL = 0.024           # ps

# Salt
MIN_CHARGED_SALT_CONC = 0.15            # M

# GCMC — from desmond/msj_constants.py _GCMC_SOLVENT_VDW_SCALE_FACTOR
GCMC_SOLVENT_VDW_SCALE_FACTOR = 0.75

# ENERGY_GROUP output — from desmond/msj_constants.py ENERGY_GROUP
ENEGRP_NAME_TEMPLATE = "$JOBNAME$[_replica$REPLICA$]_enegrp.dat"
ENEGRP_FIRST = 0.0    # ps
ENEGRP_INTERVAL = 1.2  # ps

# RBFE
RBFE_DEFAULT_LAMBDAS = 12

# ABFE lambda window counts — from desmond/fep_schedule.py:282-301
# Default (neutral) protocol
ABFE_COMPLEX_LAMBDAS = 68
ABFE_SOLVENT_LAMBDAS = 60
ABFE_RESTRAINED_COMPLEX_LAMBDAS = 80
ABFE_RESTRAINED_SOLVENT_LAMBDAS = 68
# Charged protocol
ABFE_CHARGED_COMPLEX_LAMBDAS = 108
ABFE_CHARGED_SOLVENT_LAMBDAS = 60
ABFE_CHARGED_RESTRAINED_COMPLEX_LAMBDAS = 128
ABFE_CHARGED_RESTRAINED_SOLVENT_LAMBDAS = 68

# Backbone ASL for receptor restraints — from absolute_binding.py
BACKBONE_ASL = (
    "((protein and backbone) or "
    "(nucleic_acids and nucleic_backbone)) and not a.ele H"
)
RECEPTOR_ASL = "protein or nucleic_acids"


class StrucTag:
    """Values of the s_fep_struc_tag CT property — from desmond/constants.py FEP_STRUC_TAG."""
    RECEPTOR = "receptor"
    LIGAND   = "ligand"
    COMPLEX  = "complex"
    SOLVENT  = "solvent"
    MEMBRANE = "membrane"


# Per-leg CT selection — mirrors desmond/msj_generator/keep_struc_tags.py
KEEP_STRUC_TAGS: dict[str, dict[str, list[str]]] = {
    "rbfe": {
        "complex": [
            StrucTag.RECEPTOR, StrucTag.MEMBRANE,
            StrucTag.SOLVENT,  StrucTag.LIGAND,
        ],
        "solvent": [StrucTag.LIGAND],
        "vacuum":  [StrucTag.LIGAND],
    },
    "abfe": {
        "md": [
            StrucTag.RECEPTOR, StrucTag.SOLVENT,
            StrucTag.MEMBRANE, StrucTag.LIGAND,
        ],
        "complex": [
            StrucTag.RECEPTOR, StrucTag.MEMBRANE,
            StrucTag.SOLVENT,  StrucTag.COMPLEX,
        ],
        "solvent": [StrucTag.LIGAND],
    },
}
