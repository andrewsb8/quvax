from src.qodon.optimizers.metro_optimizer import MetropolisOptimizer
import random


class REMCOptimizer(MetropolisOptimizer):
    """
    Replica Excahnge Metropolis Monte Carlo Optimizer for codon optimization

    Parameters
    -----------

    """

    def __init__(self, config):
        super().__init__(config)
        self.beta_list = self._generate_temperatures()

    def _optimize(self):
        """
        Method for codon optimization using replica exchange metropolis monte
        carlo algorithm

        """

        if not self.config.args.resume:
            self._iterate(self.initial_sequences, update_counter=False)
            members = self.initial_sequences
            sec_structs = self.sec_structs
        else:
            members = [self._convert_codons_to_ints(s) for s in self.initial_sequences]
            # not loading previous secondary structure because they are not compared
            sec_structs = ["" for i in range(self.config.args.population_size)]

        energies = self.energies
        self.accepted = 0
        self.rejected = 0
        self.randomed = 0
        self.accepted_exchanges = 0
        self.rejected_exchanges = 0
        self.exchange_attempt = 0

        for i in range(self.config.args.codon_iterations):
            if i != 0 and self.config.args.exchange_frequency % i == 0:
                self._attempt_exchanges(i, members, energies)
            self._metropolis_iteration(members, energies, sec_structs)
            self._iterate(
                members, energies, sec_structs
            )  # pass energies and ss to _iterate to avoid refolding

        if self.config.args.resume:
            self.config.log.info(
                """MC stats are only kept track of in each individual execution of
                design.py and are not aggregated from previous runs if using --resume."""
            )
        self.config.log.info("Sequence changes accepted: " + str(self.accepted))
        self.config.log.info("Sequence changes rejected: " + str(self.rejected))
        self.config.log.info("Random sequences generated: " + str(self.randomed))
        self.config.log.info(
            "System exchanges accepted: " + str(self.accepted_exchanges)
        )
        self.config.log.info(
            "System exchanges rejected: " + str(self.rejected_exchanges)
        )
        self._post_process()

    def _generate_temperatures(self):
        """
        define range of beta values in list

        """
        return [
            (
                self.config.args.beta
                + i
                * (self.config.args.beta_max - self.config.args.beta)
                / self.config.args.population_size
            )
            for i in range(self.config.args.population_size)
        ]

    def _attempt_exchanges(self, iteration, members, energies):
        """
        try to exchange systems to different temperatures

        """
        # if first digit is even, only try to exchange even indices
        if int(str(iteration)[0]) % 2 == 0:
            for i in range(1, self.config.args.population_size, 2):
                if self._check_changes(
                    self.beta_list[i - 1],
                    energies[i - 1],
                    self.beta_list[i],
                    energies[i],
                ):
                    self._accept_exchange(i, members, energies)
                else:
                    self.rejected_exchanges += 1
        else:
            for i in range(2, self.config.args.population_size, 2):
                if self._check_changes(
                    self.beta_list[i - 1],
                    energies[i - 1],
                    self.beta_list[i],
                    energies[i],
                ):
                    self._accept_exchange(i, members, energies)
                else:
                    self.rejected_exchanges += 1

    def _accept_exchange(self, i, members, energies):
        members[i-1], members[i] = members[i], members[i-1]
        energies[i-1], energies[i] = energies[i], energies[i-1]
        self.accepted_exchanges += 1
