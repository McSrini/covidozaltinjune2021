/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ca.mcmaster.covidozaltinjune2021.cplex;

import static ca.mcmaster.covidozaltinjune2021.Constants.DOUBLE_ONE;
import static ca.mcmaster.covidozaltinjune2021.Constants.ONE;
import static ca.mcmaster.covidozaltinjune2021.Constants.ZERO;
import ca.mcmaster.covidozaltinjune2021.Parameters;
import static ca.mcmaster.covidozaltinjune2021.Parameters.K_I_J;
import static ca.mcmaster.covidozaltinjune2021.Parameters.NUMBER_OF_SUBGROUPS_J;
import static ca.mcmaster.covidozaltinjune2021.cplex.Variable_Creator.createContinuousVars;
import ilog.concert.IloException;
import ilog.concert.IloNumExpr;
import ilog.concert.IloNumVar;
import ilog.cplex.IloCplex;

/**
 *
 * @author tamvadss
 */
public class LP_Model_For_Bound_Initialization {
    
    private  IloCplex  cplex;
    private IloNumVar[] v_Vars = null;
    
    public   double getBound (int ii, boolean isUpper) throws IloException{
       cplex= new IloCplex() ; 
       
       v_Vars = createContinuousVars(cplex, "v", NUMBER_OF_SUBGROUPS_J, DOUBLE_ONE) ;
       cplex.add(v_Vars);
       
       addObjective(ii, isUpper ) ;
       addConstraint_6C();
       
       IloNumExpr constraint = cplex.numExpr();
       for (int jj = ZERO; jj < NUMBER_OF_SUBGROUPS_J; jj ++){
            constraint =cplex .sum (constraint, cplex.prod  (K_I_J.get(ii).get(jj),v_Vars[jj] ) ) ;
       }
       constraint = cplex.prod (constraint, DOUBLE_ONE - Parameters.VACCINE_EFFICACY_PSI );
       cplex.addLe (constraint, v_Vars[ii]) ;
       
       cplex.exportModel("LP_for_bounds.lp");
       
       cplex.solve ();
       
       return cplex.getObjValue();
    }
    
    private void addObjective (int ii, boolean isUpperBound) throws IloException{
        IloNumExpr objective = cplex.numExpr();
        for (int jj = ZERO; jj < NUMBER_OF_SUBGROUPS_J; jj ++){
            objective =cplex .sum (objective, cplex.prod  (K_I_J.get(ii).get(jj),v_Vars[jj] ) ) ;
        }
        if (isUpperBound){
            cplex.addMaximize (objective );            
        }else {
            cplex.addMinimize (objective );
        }
    }
    
    private   void addConstraint_6C( ) throws IloException{
        IloNumExpr numericExpr  = cplex.numExpr();
        for (IloNumVar var : v_Vars){
            numericExpr= cplex.sum ( numericExpr , var);
        }
        cplex.addEq ( numericExpr, ONE);
    }
}
