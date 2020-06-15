# MizerEvoPlugin
This project adds the possibility to simulate evolutionary processes with the Mizer model

At the moment the code is available but needs to be linked with the base mizer package.
The code uses function from `mizer` and doesn't work if `mizer` is simply loaded with `library()`. `mizer`'s functions needs to be available in the environment.

What the code can do:

- Insert a predetermined number of new species (`mutation` parameter) during a mizer projection (`evoProject()`)

- New species are copy of existing ones albeit with a change to their maturation size

- Species reaching an abundance lower than $10^{-30}$ are set to 0 abundance with `extinctionRDD()`

- Each projections containing a new species is independent of the others. They are saved in a temporary folder before being binded at the end of the simulation