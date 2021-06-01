/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ca.mcmaster.covidozaltinjune2021.cplex;

import static ca.mcmaster.covidozaltinjune2021.Constants.*; 
import ca.mcmaster.covidozaltinjune2021.Parameters;
import static ca.mcmaster.covidozaltinjune2021.Parameters.*;
import static ca.mcmaster.covidozaltinjune2021.cplex.Variable_Creator.*;
import ilog.concert.IloException;
import ilog.concert.IloNumExpr;
import ilog.concert.IloNumVar;

/**
 *
 * @author tamvadss
 */
public class Lower_ModelBuilder extends Base_ModelBuilder {
    
    protected IloNumVar[] s_Vars=null;
    protected IloNumVar[] h_Vars=null;
    
    public Lower_ModelBuilder() throws IloException{
        super();
        
        s_Vars = createContinuousVars(cplex, "s", NUMBER_OF_SUBGROUPS_J) ;
        h_Vars = createContinuousVars(cplex, "h", NUMBER_OF_SUBGROUPS_J, DOUBLE_BILLION) ;
    }

    @Override
    protected void buildModel() throws IloException {
        //        
        for (int ii = ZERO; ii < NUMBER_OF_SUBGROUPS_J; ii ++){
            IloNumVar[] xi =  createBinaryVars(cplex, "x_"+ ii+"_", ONE+ Parameters.BIN_VARS_FOR_DISCRETIZATION_P) ;
            x_VarList.add (xi);
            cplex.add(xi);
            addConstraint_10( ii );
        }
        
        for (int ii = ZERO; ii < NUMBER_OF_SUBGROUPS_J; ii ++){
            IloNumVar[] wi =  createContinuousVars(cplex, "w_"+ ii+"_",  ONE + Parameters.BIN_VARS_FOR_DISCRETIZATION_P, DOUBLE_ONE) ;
            w_VarList.add (wi);
            cplex.add(wi);             
        }
        for (int ii = ZERO; ii < NUMBER_OF_SUBGROUPS_J; ii ++){
            for (int el = ONE; el < ONE+ BIN_VARS_FOR_DISCRETIZATION_P; el ++){                
                cplex.addLe(w_VarList.get(ii)[el],  x_VarList.get(ii)[el]);
                cplex.addLe( w_VarList.get(ii)[el],  f_Vars [ii] );
            }
        }
             
        for (int ii = ZERO; ii < NUMBER_OF_SUBGROUPS_J; ii ++){
            addConstraint_16 (  ii );
            addConstraint_12a (  ii );
        }
       
        cplex.exportModel( "LB.lp");
                
    }
    
    
    private  void addConstraint_16 (int ii ) throws IloException{
        
        cplex.addLe( h_Vars[ii], s_Vars[ii] );
        cplex.addLe( h_Vars[ii], cplex.prod( Math.pow( DOUBLE_TWO, -Parameters.BIN_VARS_FOR_DISCRETIZATION_P ),  f_Vars[ii] ));
        
        IloNumExpr numericExpr   = cplex.numExpr();
        numericExpr= cplex.sum (numericExpr, s_Vars[ii] );
        numericExpr= cplex.sum (numericExpr, -Math.pow(DOUBLE_TWO, -BIN_VARS_FOR_DISCRETIZATION_P  ) );
        IloNumExpr numericExpr_prod   = cplex.numExpr();
        numericExpr_prod = cplex.sum (numericExpr_prod,f_Vars[ii] );
        numericExpr_prod= cplex.prod(numericExpr_prod, Math.pow(DOUBLE_TWO, -BIN_VARS_FOR_DISCRETIZATION_P  )  );
        numericExpr= cplex.sum (numericExpr, numericExpr_prod );
        cplex.addLe(numericExpr,  h_Vars[ii]);
        
    }
    
    private   void addConstraint_10(int ii ) throws IloException{
        IloNumExpr numericExpr_LHS  = cplex.numExpr();
        
        for (int jj = ZERO; jj <NUMBER_OF_SUBGROUPS_J ; jj ++){
            numericExpr_LHS = cplex.sum (numericExpr_LHS,cplex.prod(v_Vars[jj] , Parameters.K_I_J.get(ii).get(jj) ) ) ;
        }
        
        numericExpr_LHS = cplex.sum (numericExpr_LHS,- LL_I.get(ii) ) ;
        numericExpr_LHS = cplex.sum(numericExpr_LHS,cplex.prod(-ONE,  s_Vars[ii])) ;
        
        numericExpr_LHS = cplex.prod (numericExpr_LHS,DOUBLE_ONE/( UL_I.get(ii)- LL_I.get(ii))) ;
        
        IloNumExpr numericExpr_RHS  = cplex.numExpr();
        for (int el = ZERO; el< ONE+ Parameters.BIN_VARS_FOR_DISCRETIZATION_P; el++){
            numericExpr_RHS = cplex.sum (numericExpr_RHS, cplex.prod (Math.pow( DOUBLE_TWO, -el), x_VarList.get(ii)[el])) ;
        }
        
        
        
        cplex.addEq(numericExpr_LHS, numericExpr_RHS);
    }
    
    
    private   void addConstraint_12a(int ii ) throws IloException{
        IloNumExpr numericExpr_ONE  = cplex.numExpr();
        IloNumExpr numericExpr_TWO  = cplex.numExpr();
        IloNumExpr numericExpr_THREE = cplex.numExpr();
        IloNumExpr numericExpr_LHS = cplex.numExpr();
        
        numericExpr_ONE= cplex.prod (f_Vars[ii], Parameters.VACCINE_EFFICACY_PSI );
        numericExpr_ONE = cplex.prod(numericExpr_ONE, -ONE);
        numericExpr_ONE = cplex.sum( ONE,numericExpr_ONE );
        numericExpr_ONE = cplex.prod(numericExpr_ONE,LL_I.get(ii) );
        
        for (int el = ZERO; el< ONE+ Parameters.BIN_VARS_FOR_DISCRETIZATION_P; el++){
            numericExpr_TWO = cplex.sum (numericExpr_TWO, cplex.prod (Math.pow( DOUBLE_TWO, -el), x_VarList.get(ii)[el])) ;
        }
        numericExpr_TWO= cplex.prod(numericExpr_TWO, UL_I.get(ii) - LL_I.get(ii));
        
        
        
        for (int el = ZERO; el< ONE+ Parameters.BIN_VARS_FOR_DISCRETIZATION_P; el++){
            numericExpr_THREE = cplex.sum (numericExpr_THREE, cplex.prod (Math.pow( DOUBLE_TWO, -el), w_VarList.get(ii)[el])) ;
        }
        numericExpr_THREE= cplex.prod(numericExpr_THREE, UL_I.get(ii) - LL_I.get(ii));
        numericExpr_THREE= cplex.prod(numericExpr_THREE,  - Parameters.VACCINE_EFFICACY_PSI);        
        
        
        numericExpr_LHS = cplex.sum (numericExpr_LHS, numericExpr_ONE);
        numericExpr_LHS = cplex.sum (numericExpr_LHS, numericExpr_TWO);
        numericExpr_LHS = cplex.sum (numericExpr_LHS, numericExpr_THREE);
        
        
        IloNumExpr numericExpr_LHS_addendum = cplex.numExpr();
        ///add si -phi*hi
        numericExpr_LHS_addendum = cplex.sum (numericExpr_LHS_addendum, s_Vars[ii]) ;
        numericExpr_LHS_addendum = cplex.sum (numericExpr_LHS_addendum, 
                cplex.prod(-Parameters.VACCINE_EFFICACY_PSI,h_Vars[ii] )) ;
        
        numericExpr_LHS= cplex.sum (numericExpr_LHS, numericExpr_LHS_addendum );
        cplex.addLe(  numericExpr_LHS, v_Vars[ii]);
        
    }
    
}
