from src.qodon.optimizer import CodonOptimizer

class MockOptimizer(CodonOptimizer):
    def __init__(self, config):
        super().__init__(config)

    def _optimize(self):
        pass
