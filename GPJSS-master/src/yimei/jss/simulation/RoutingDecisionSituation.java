package yimei.jss.simulation;

import yimei.jss.jobshop.OperationOption;
import yimei.jss.simulation.state.SystemState;

import java.util.ArrayList;
import java.util.List;

/**
 * Created by dyska on 7/09/17.
 */
public class RoutingDecisionSituation extends DecisionSituation {

    private List<OperationOption> queue;
    //private WorkCenter workCenter;
    private SystemState systemState;

    public RoutingDecisionSituation(List<OperationOption> operationOptions,
                                    //WorkCenter workCenter,
                                    SystemState systemState) {
        this.queue = operationOptions;
        //this.workCenter = workCenter;
        this.systemState = systemState;
    }

    public List<OperationOption> getQueue() {
        return queue;
    }
    /*public WorkCenter getWorkCenter() {
        return workCenter;
    }*/
    public SystemState getSystemState() {
        return systemState;
    }


    public RoutingDecisionSituation clone() {

        //Yi
        List<OperationOption> clonedQ = new ArrayList<>();
        for (OperationOption op : queue) {
            clonedQ.add(op.clone());
        }
        SystemState clonedState = systemState.clone();

      /*  List<OperationOption> clonedQ = new ArrayList<>(queue);
        WorkCenter clonedWC = workCenter.clone();
        SystemState clonedState = systemState.clone();
        return new RoutingDecisionSituation(clonedQ, clonedWC, clonedState);*/

        return new RoutingDecisionSituation(clonedQ, clonedState);
    }
}
