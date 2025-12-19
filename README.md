# Synaptic potentiation dependence on spike variability

## INTRODUCTION
Understanding how the brain modifies synaptic connections and how firing patterns influence this process remains a major challenge in neuroscience. Aiming to clarify this relationship, we consider all-to-all neuronal networks with spike timing-dependent plasticity (STDP) where neurons are connected through excitatory chemical synapses with weak initial synaptic weights. We analyze how spike variability within a phase-synchronous pattern affects synaptic potentiation between neurons. Considering different methodologies, we find that, depending on the variability of spike synchronization and firing frequency, the potentiation of neuronal connections generates predominant unidirectional or bidirectional topologies. In addition, we identify four types of triad structures that are induced in the network. Particularly, for a certain level of variability in phase synchronization, a non-trivial optimization of the mean potentiation per spike is observed. In these cases, the potentiation occurs at a higher rate due to the preferential formation of unidirectional connections. Overall, our results deepen the knowledge of how phase firing patterns drive the synaptic changes in neuronal networks in the presence of STDP.

## COMPILE AND RUN
Compile and run using:
\`\`\`bash
bash compilar_rodar.sh
\`\`\`

Join the parallelized files using:
\`\`\`bash
join_files.bash
\`\`\`

## PLOTTING
To plot Fig. 4 use:
\`\`\`bash
gnuplot fig4.plt
\`\`\`

Simulation and data of Fig.5 is avaliable in file:
\`\`\`
Fig5
\`\`\`

To plot Fig. 6 use:
\`\`\`bash
python3 fig5.py
\`\`\`

## FINAL NOTE
Additional questions and data can be addressed to:
\`\`\`
protachevicz@gmail.com
\`\`\`
