package yimei.jss.algorithm.multitreeModelSurrogateSample;

import ec.EvolutionState;
import ec.Individual;
import ec.gp.GPIndividual;
import ec.gp.GPNode;
import ec.simple.SimpleProblemForm;
import ec.simple.SimpleStatistics;
import yimei.jss.gp.GPRuleEvolutionState;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import static yimei.jss.algorithm.multitreeModelSurrogateSample.PhenotypicGPRuleEvolutionState.indexEvaluated;
import static yimei.jss.gp.GPRun.out_dir;

public class SimpleStatisticsSaveRulesizePrintInds extends SimpleStatistics {

    //get seed
    protected long jobSeed;

    //fzhang 25.6.2018 in order to save the rulesize in each generation
    List<Long> aveSeqRulesizeTree0 = new ArrayList<>();
    List<Long> aveRouRulesizeTree1 = new ArrayList<>();

    /**
     * GENERATIONAL: Called immediately before evaluation occurs.
     */
    public void preEvaluationStatistics(final EvolutionState state) {
        for (int x = 0; x < children.length; x++)
            children[x].preEvaluationStatistics(state);

        //fzhang 17.6.2018  get the seed value
//    	Parameter p;
//  		// Get the job seed.
//  		p = new Parameter("seed").push(""+0);
//        jobSeed = state.parameters.getLongWithDefault(p, null, 0);

        // fzhang 15.6.2018 1. save the individual size in population
        // 2. calculate the average size of individuals in population
        // check the average size of sequencing and routing rules in population
        //fzhang 15.6.2018  in order to check the average size of sequencing and routing rules in population
        int SeqSizeTree0 = 0;
        int RouSizeTree1 = 0;
        //int indSizePop = 0; // in order to check whether SeqSizePop1 and RouSizePop2 are calculated correctly
        // should be the sum of SeqSizePop1 and RouSizePop2
        long aveSeqSizeTree0 = 0;
        long aveRouSizeTree1 = 0;
        for (int ind = 0; ind < state.population.subpops[0].individuals.length; ind++) {
            GPIndividual indi = (GPIndividual) state.population.subpops[0].individuals[ind];
            SeqSizeTree0 += indi.trees[0].child.numNodes(GPNode.NODESEARCH_ALL);
            RouSizeTree1 += indi.trees[1].child.numNodes(GPNode.NODESEARCH_ALL);
        }
        aveSeqSizeTree0 = SeqSizeTree0 / state.population.subpops[0].individuals.length;
        aveRouSizeTree1 = RouSizeTree1 / state.population.subpops[0].individuals.length;

        aveSeqRulesizeTree0.add(aveSeqSizeTree0);
        aveRouRulesizeTree1.add(aveRouSizeTree1);

//                if (state.generation == state.numGenerations - 1) {
        if (((GPRuleEvolutionState)state).totalTime >= ((GPRuleEvolutionState)state).endTime) {
            //fzhang  15.6.2018  save the size of rules in each generation
            File rulesizeFile = new File(out_dir + "/job." + GPRuleEvolutionState.jobSeed + ".aveGenRulesize.csv"); // jobSeed = 0

            try {
                BufferedWriter writer = new BufferedWriter(new FileWriter(rulesizeFile));
                writer.write("Gen,aveSeqRuleSize,aveRouRuleSize,avePairSize");
                writer.newLine();
                for (int gen = 0; gen < aveSeqRulesizeTree0.size(); gen++) {
                    writer.write(gen + "," + aveSeqRulesizeTree0.get(gen) + "," + aveRouRulesizeTree1.get(gen) + "," +
                            (aveSeqRulesizeTree0.get(gen) + aveRouRulesizeTree1.get(gen)) / 2);
                    writer.newLine();
                }
                writer.close();
            } catch (IOException e) {
                e.printStackTrace();
            }

			/*System.out.println(SeqSizeTree0);
			System.out.println(RouSizeTree1);
			System.out.println(aveSeqRulesizeTree0.get(state.generation));
			System.out.println(aveRouRulesizeTree1.get(state.generation));*/

            //fzhang 15.6.2018 in order to check whether SeqSizePop1 and RouSizePop2 are calculated correctly (YES)
	 	/*	for (int pop = 0; pop < state.population.subpops.length; pop++) {
	 			for (int ind = 0; ind < state.population.subpops[pop].individuals.length; ind++) {
	 				indSizePop += state.population.subpops[pop].individuals[ind].size();
	 			}
	 		}
	 		System.out.println(indSizePop);*/
        }
    }

    boolean warned = false;

    public void postEvaluationStatistics(final EvolutionState state) {

        // for now we just print the best fitness per subpopulation.
        Individual[] best_i = new Individual[state.population.subpops.length];  // quiets compiler complaints
        for (int x = 0; x < state.population.subpops.length; x++) {
            best_i[x] = state.population.subpops[x].individuals[indexEvaluated.get(0)];
            for (int a = 1; a < indexEvaluated.size(); a++) {
                int y = indexEvaluated.get(a);
                if (state.population.subpops[x].individuals[y] == null) {
                    if (!warned) {
                        state.output.warnOnce("Null individuals found in subpopulation");
                        warned = true;  // we do this rather than relying on warnOnce because it is much faster in a tight loop
                    }
                } else if (best_i[x] == null || state.population.subpops[x].individuals[y].fitness.betterThan(best_i[x].fitness))
                    best_i[x] = state.population.subpops[x].individuals[y];
                if (best_i[x] == null) {
                    if (!warned) {
                        state.output.warnOnce("Null individuals found in subpopulation");
                        warned = true;  // we do this rather than relying on warnOnce because it is much faster in a tight loop
                    }
                }
            }

            // now test to see if it's the new best_of_run
            if (best_of_run[x] == null || best_i[x].fitness.betterThan(best_of_run[x].fitness))
                best_of_run[x] = (Individual) (best_i[x].clone());
        }

        // print the best-of-generation individual
        if (doGeneration) state.output.println("\nGeneration: " + state.generation, statisticslog);
        if (doGeneration) state.output.println("Best Individual:", statisticslog);
        for (int x = 0; x < state.population.subpops.length; x++) {
            if (doGeneration) state.output.println("Subpopulation " + x + ":", statisticslog);
            if (doGeneration) best_i[x].printIndividualForHumans(state, statisticslog);
            if (doMessage && !silentPrint) state.output.message("Subpop " + x + " best fitness of generation" +
                    (best_i[x].evaluated ? " " : " (evaluated flag not set): ") +
                    best_i[x].fitness.fitnessToStringForHumans());

            // describe the winner if there is a description
            if (doGeneration && doPerGenerationDescription) {
                if (state.evaluator.p_problem instanceof SimpleProblemForm)
                    ((SimpleProblemForm) (state.evaluator.p_problem.clone())).describe(state, best_i[x], x, 0, statisticslog);
            }
        }

        indexEvaluated.clear();
    }

}
