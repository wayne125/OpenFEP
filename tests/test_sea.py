from openfep.msj.sea import SEABlock, render, sea_val, render_msj


def test_sea_val_string():
    assert sea_val("OPLS4") == '"OPLS4"'


def test_sea_val_bool_true():
    assert sea_val(True) == "true"


def test_sea_val_bool_false():
    assert sea_val(False) == "false"


def test_sea_val_int():
    assert sea_val(42) == "42"


def test_sea_val_float():
    assert sea_val(3.14) == "3.14"


def test_sea_val_list_of_numbers():
    assert sea_val([0.004, 0.004, 0.008]) == "[0.004 0.004 0.008]"


def test_simple_block_renders():
    block = SEABlock("assign_forcefield", ("forcefield", "OPLS4"))
    text = render(block)
    assert text.startswith("assign_forcefield {")
    assert 'forcefield = "OPLS4"' in text
    assert text.endswith("}")


def test_block_bool_field_lowercase():
    block = SEABlock("assign_forcefield", ("hydrogen_mass_repartition", True))
    text = render(block)
    assert "hydrogen_mass_repartition = true" in text


def test_nested_block_renders():
    inner = SEABlock("simulate", ("fep.lambda", "default:12"))
    outer = SEABlock("set_family", ("simulate", inner))
    text = render(outer)
    assert "set_family {" in text
    assert "simulate {" in text
    assert 'fep.lambda = "default:12"' in text


def test_render_msj_joins_with_blank_line():
    a = SEABlock("task", ("task", "generic"))
    b = SEABlock("simulate", ("time", 120.0))
    text = render_msj([a, b])
    assert "task {" in text
    assert "simulate {" in text
    assert "\n\n" in text   # blank line between stages


def test_indentation():
    block = SEABlock("simulate", ("time", 120.0))
    text = render(block, indent=1)
    assert text.startswith("  simulate {")   # 2 spaces per indent level
    assert "    time = 120.0" in text         # child indented by 4
