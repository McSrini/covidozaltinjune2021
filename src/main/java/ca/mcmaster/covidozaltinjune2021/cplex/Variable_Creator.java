/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ca.mcmaster.covidozaltinjune2021.cplex;

import static ca.mcmaster.covidozaltinjune2021.Constants.DOUBLE_TWO;
import static ca.mcmaster.covidozaltinjune2021.Constants.ZERO;
import ca.mcmaster.covidozaltinjune2021.Parameters;
import ilog.concert.IloException;
import ilog.concert.IloNumVar;
import ilog.cplex.IloCplex;

/**
 *
 * @author tamvadss
 */
public class Variable_Creator {
    
        
    public static IloNumVar[] createContinuousVars (IloCplex cplex, String prefix, int count) throws IloException{
        String[] xName = new String[count];
        for (int index = ZERO; index < count; index ++){
            xName[index] = prefix+index;
        }
        IloNumVar[] createdVars=  cplex.numVarArray( count, ZERO, Math.pow (DOUBLE_TWO, -Parameters.BIN_VARS_FOR_DISCRETIZATION_P), xName);
        return createdVars;
    }
    
    public static IloNumVar[] createContinuousVars (IloCplex cplex, String prefix, int count, Double UPPER_BOUND) throws IloException{
        String[] xName = new String[count];
        for (int index = ZERO; index < count; index ++){
            xName[index] = prefix+index;
        }
        IloNumVar[] createdVars=  cplex.numVarArray( count, ZERO,  UPPER_BOUND , xName);
        return createdVars;
    }
    
     
    public static IloNumVar[] createBinaryVars (IloCplex cplex, String prefix, int count) throws IloException{
        String[] xName = new String[count];
        for (int index = ZERO; index < count; index ++){
            xName[index] = prefix+index;
        }
        IloNumVar[] createdVars=  cplex.boolVarArray (count,    xName) ;
        
        return createdVars;
    }
     
    
}
