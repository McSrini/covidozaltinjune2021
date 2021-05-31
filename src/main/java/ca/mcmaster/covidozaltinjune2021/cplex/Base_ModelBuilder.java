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
import static ca.mcmaster.covidozaltinjune2021.Parameters.NUMBER_OF_SUBGROUPS_J;
import static ca.mcmaster.covidozaltinjune2021.cplex.Variable_Creator.*;
import ilog.concert.IloException;
import ilog.concert.IloNumExpr;
import ilog.concert.IloNumVar;
import ilog.cplex.IloCplex;
import static ilog.cplex.IloCplex.Status.Optimal;
import java.util.ArrayList;
import java.util.List;

/**
 *
 * @author tamvadss
 */
public abstract class Base_ModelBuilder {
    
    protected  IloCplex cplex  ;
    protected IloNumVar[] v_Vars=null;
    protected IloNumVar[] f_Vars=null;
       
    protected     List<Double> UL_I = new ArrayList<Double>() ;
    protected      List<Double> LL_I = new ArrayList<Double>() ;
    
    protected List< IloNumVar[]>  x_VarList = new ArrayList < IloNumVar[]> ();    
    protected List< IloNumVar[]>  w_VarList = new ArrayList < IloNumVar[]> ();
    
      
    public Base_ModelBuilder () throws IloException{
        
        cplex= new IloCplex() ;
        
        //init upper and lower bounds
        UL_I.addAll( Parameters.UL_I);
        LL_I.addAll( Parameters.LL_I);
        
            
        v_Vars = createContinuousVars(cplex, "v", NUMBER_OF_SUBGROUPS_J, DOUBLE_ONE) ;
        f_Vars = createContinuousVars(cplex, "f", NUMBER_OF_SUBGROUPS_J, DOUBLE_ONE) ;
        
        cplex.add(v_Vars);
        cplex.add(f_Vars);
         
        addConstraint_6C();
        
        addMinObjective();
        
        cplex.exportModel("base.lp" );
        
    }
     
    
    public void solve () throws IloException{                
          
        buildModel();

        //solve both upper and lower cplex
        cplex.solve();
        System.out.println ("cplex.getStatus() "+ cplex.getStatus()) ;
        if (cplex.getStatus().equals(Optimal)){
            System.out.println("\n\nBound is : "+cplex.getObjValue() + "\n\n" );
        }
        
        
                
    }
    
    protected void addMinObjective () throws IloException{
        IloNumExpr objective = cplex.numExpr();
        for (int ii = ZERO; ii < NUMBER_OF_SUBGROUPS_J; ii ++){
            objective =cplex .sum (objective, cplex.prod (f_Vars[ii], Parameters.SUBGROUP_SIZES.get(ii)) ) ;
        }
        cplex.addMinimize (objective );
    }
    
    protected   void addConstraint_6C( ) throws IloException{
        IloNumExpr numericExpr  = cplex.numExpr();
        for (IloNumVar var : v_Vars){
            numericExpr= cplex.sum ( numericExpr , var);
        }
        cplex.addEq ( numericExpr, ONE);
    }
   
     
    protected abstract    void buildModel ( ) throws IloException;
    
}
