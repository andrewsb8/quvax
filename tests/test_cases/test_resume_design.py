import shutil
from src.params.design_parser import DesignParser


def test_resume(caplog):
    """
    Test to verify successful execution of resuming optimization

    """
    from src.qodon.optimizers.classical_ga import GeneticAlgorithm

    #copy database so test does not add information to test file
    #it will be deleted after test is completed
    shutil.copy("tests/test_files/test_design/quvax.db", "quvax.db")
    testargs = [
        "-i",
        "quvax.db",
        "--resume",
    ]
    config = DesignParser._resume(testargs)
    GeneticAlgorithm(config)
    log_entry = (
        "src.params.design_parser",
        20,  # 40 indicates error, 30 indicates WARNING, 20 indicates INFO
        "Finished parsing optimized sequences.",
    )
    assert log_entry in caplog.record_tuples
