# Scalable Spike-and-Slab MCMC
# Biswas, Mackey and Meng, ICML 2022.

from enum import Enum
import pandas as pd
from typing import Dict, Optional
from mcmc_functions import spike_slab_linear, spike_slab_logistic, spike_slab_probit

# Model type constants
class ModelType(Enum):
    LINEAR = "linear"
    LOGISTIC = "logistic"
    PROBIT = "probit"


class SpikeSlabPrior:
    def __init__(
        self,
        model_type: ModelType,
        q: float,
        tau0: float,
        tau1: float,
        a0: Optional[float] = None,
        b0: Optional[float] = None,
    ):
        self.model_type = model_type
        self.q = q
        self.tau0 = tau0
        self.tau1 = tau1
        self.a0 = a0
        self.b0 = b0

    def get_params(self):
        if self.model_type == ModelType.LINEAR:
            if (self.a0 is None) or (self.b0 is None):
                ValueError(
                    f"{ModelType.LINEAR} models cannot have a0 and b0 parameters as None"
                )
            return {
                "q": self.q,
                "tau0": self.tau0,
                "tau1": self.tau1,
                "a0": self.a0,
                "b0": self.b0,
            }
        elif self.model_type in (ModelType.LOGISTIC, ModelType.PROBIT):
            if (self.a0 is not None) or (self.b0 is not None):
                ValueError(
                    f"{ModelType.LOGISTIC.value} or {ModelType.PROBIT.value} models must have a0 and b0 parameters as None"
                )
            return {"q": self.q, "tau0": self.tau0, "tau1": self.tau1}


class ScaleSpikeSlab:
    def __init__(
        self,
        endog: pd.Series,
        exog: pd.DataFrame,
        model: ModelType,
        prior: SpikeSlabPrior,
    ):
        self.endog = endog
        self.exog = exog
        self.model = model
        self.prior = prior
        self.prior = prior

    def run_mcmc(
        self,
        chain_length: int,
        burnin: int,
        store_trajectory: bool = False,
        rinit=None,
        verbose=False,
    ):
        if self.model == ModelType.LINEAR:
            chain_output = spike_slab_linear(
                y=self.endog,
                X=self.exog,
                rinit=rinit,
                verbose=verbose,
                chain_length=chain_length,
                burnin=burnin,
                store=store_trajectory,
                **self.prior.get_params(),
            )
        elif self.model == ModelType.LOGISTIC:
            chain_output = spike_slab_logistic(
                y=self.endog,
                X=self.exog,
                rinit=rinit,
                verbose=verbose,
                chain_length=chain_length,
                burnin=burnin,
                store=store_trajectory,
                **self.prior.get_params(),
            )
        elif self.model == ModelType.PROBIT:
            chain_output = spike_slab_probit(
                y=self.endog,
                X=self.exog,
                rinit=rinit,
                verbose=verbose,
                chain_length=chain_length,
                burnin=burnin,
                store=store_trajectory,
                **self.prior.get_params(),
            )
        return chain_output
