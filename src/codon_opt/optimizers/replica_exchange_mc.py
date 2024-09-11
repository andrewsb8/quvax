from src.codon_opt.optimizers.metro_optimizer import MetropolisOptimizer
import random


class REMCOptimizer(MetropolisOptimizer):
    """
    Replica Excahnge Metropolis Monte Carlo Optimizer for codon optimization.
    This is a direct extension of the Metropolis Monte Carlo Optimizer where
    each member of the population is exposed to a temperature within a range
    set by the user. Periodically, at a frequency set by the user, the sequences
    are checked to see if they can be swapped between temperatures using a similar
    Boltzmann criterion governing the Metropolis algorithm (except with different
    beta (1/kT) values).

    Parameters
    -----------
    self.beta_list : list float
        List of temperatures for each sequence in the population which has range
        self.config.args.beta to self.config.args.beta_max and the interval is
        the difference of those values divided by population size

    """

    def __init__(self, config):
        super().__init__(config)
        self.beta_list = self._generate_temperatures()
        self.accepted_exchanges = 0
        self.rejected_exchanges = 0
        self.exchange_attempt = 0

    def _optimize(self):
        """
        Method for codon optimization using replica exchange metropolis monte
        carlo algorithm

        """

        for i in range(self.config.args.codon_iterations):
            if (
                i != 0
                and i != self.config.args.codon_iterations - 1
                and self.config.args.exchange_frequency % i == 0
            ):
                self.exchange_attempt += 1
                self._attempt_exchanges(
                    self.exchange_attempt, self.members, self.energies
                )
            self._metropolis_iteration(self.members, self.energies, self.sec_structs)
            self._iterate(
                fold_sequences=False
            )  # do not calculate folding energy and secondary structures

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
        Define range of beta (1/kT) values in a list. The indices of this list
        will correspond with sequences and energies of the same indices in
        ```members``` and ```energies```

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
        Try to exchange systems to different temperatures. If first digit of
        iteration is even, only try to exchange even values of i with i+1. For
        odd, try to exchange odd values of i with i+1. This strategy is inspired
        by GROMACS (see here:
        https://manual.gromacs.org/current/reference-manual/algorithms/replica-exchange.html)

        """
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
        """
        If systems can be exchanged, then swap the sequences and energies in the
        lists so they will be associated with a new temperature in self.beta_list

        """
        members[i - 1], members[i] = members[i], members[i - 1]
        energies[i - 1], energies[i] = energies[i], energies[i - 1]
        self.accepted_exchanges += 1
