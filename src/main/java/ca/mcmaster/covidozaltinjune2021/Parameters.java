/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ca.mcmaster.covidozaltinjune2021;
 
import static ca.mcmaster.covidozaltinjune2021.Constants.*;
import ca.mcmaster.covidozaltinjune2021.cplex.LP_Model_For_Bound_Initialization;
import ilog.concert.IloException;
import static java.lang.System.exit;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author tamvadss
 */
public class Parameters {
    
    public static final int NUMBER_OF_SUBGROUPS_J = 5;
    public static final List<Integer> SUBGROUP_SIZES = Arrays.asList( 77 , 241 , 375 , 204 , 103);
    public static final double VACCINE_EFFICACY_PSI = 0.5;
    public static final int BIN_VARS_FOR_DISCRETIZATION_P = 5;
    
    
    public static final List<Double> GAMMA = Arrays.asList(1.0, 1.0, 1.0, 0.85, 0.75  );
    public static final List<Double> DELTA = Arrays.asList( 0.835,0.835,0.835,0.835,0.835 );     
    public static final List<Double> OMEGA = Arrays.asList(0.526,0.526,0.526,0.526,0.526  );    
    public static final List<Double> ALPHA = Arrays.asList(  0.435 , 0.454 , 0.327 , 0.327 , 0.327);
    public static final List<Double> LAMBDA = Arrays.asList( 0.0,0.0,0.0,0.0,0.0 );    
    public static final List<Double> proportionOfIsolatedIndividuals = Arrays.asList( 0.1,0.1,0.1,0.1,0.1 );
    
    public static  List<List<Double>> BETA_I_J = new ArrayList<List<Double>>() ;  
   
    
    ///parameters KIJ , ULi, and LL_i are computed using the other parameters
    public static final List<List<Double>> K_I_J = new ArrayList<List<Double>>() ;    
    public static final List<Double> UL_I = new ArrayList<Double>() ;
    public static final List<Double> LL_I = new ArrayList<Double>() ;
    
    static {
        if (SUBGROUP_SIZES.size()!= NUMBER_OF_SUBGROUPS_J){
            System.err.println ("NUMBER_OF_SUBGROUPS_J not equals lenghth of SUBGROUP_SIZES" );
            exit (ONE);
        }
        
        BETA_I_J.add( Arrays.asList(  0.305, 0.132 ,0.205 ,0.099 ,0.041 ));
        BETA_I_J.add( Arrays.asList(  0.032 ,0.923 ,0.158 ,0.074 ,0.028 ));
        BETA_I_J.add( Arrays.asList(  0.042, 0.132, 0.183 ,0.099 ,0.041 ));
        BETA_I_J.add( Arrays.asList(  0.032, 0.101, 0.158, 0.067 ,0.029 ));
        BETA_I_J.add( Arrays.asList(  0.032, 0.101, 0.158, 0.074 ,0.032 ));
    
        
        for (int ii = ZERO; ii< NUMBER_OF_SUBGROUPS_J; ii++){
            K_I_J.add (new ArrayList<Double>()) ;              
        }
        
        try {
            init_KIJ ();
        } catch (IloException ex) {
            System.err.println("Unable to init KIJ matrix "+ ex );
            exit(ONE);
        }
    }
    
    private static void init_KIJ () throws IloException{
        //
        for (int ii = ZERO; ii< NUMBER_OF_SUBGROUPS_J; ii++)        {
            List<Double> thisList = K_I_J.get(ii);
            for (int jj = ZERO; jj< NUMBER_OF_SUBGROUPS_J; jj++)     {
                Double firstTerm=null ;
                firstTerm = GAMMA.get(ii)* LAMBDA.get(jj)/OMEGA.get(jj);
                firstTerm = firstTerm * BETA_I_J.get(ii).get(jj);
                
                Double secondTerm=null ;
                secondTerm = GAMMA.get(ii)* DELTA.get(jj)/ALPHA.get(jj);
                secondTerm = secondTerm * BETA_I_J.get(ii).get(jj);
                secondTerm= secondTerm * (ONE -proportionOfIsolatedIndividuals.get(jj) ) ;
                
                thisList.add (firstTerm+secondTerm) ;
            }
        }
        
        for (int ii = ZERO; ii< NUMBER_OF_SUBGROUPS_J; ii++)        {
            for (int jj = ZERO; jj< NUMBER_OF_SUBGROUPS_J; jj++)        {
                System.out.print(K_I_J.get(ii).get(jj) + ", ") ;
            }
            System.out.println();
        }
        
        for (int ii = ZERO; ii< NUMBER_OF_SUBGROUPS_J; ii++)  {
            LP_Model_For_Bound_Initialization lpModel = new LP_Model_For_Bound_Initialization ();
            UL_I.add (lpModel. getBound (  ii, true) ) ;
            lpModel = new LP_Model_For_Bound_Initialization ();
            LL_I.add ( lpModel. getBound (  ii, false)  );
        }
            
        System.out.println("Ui matrix:") ;
        for (Double uli : UL_I){
            System.out.print(uli+ ", ") ;
        }
        System.out.println("\nLi matrix:") ;
        for (Double lli : LL_I){
            System.out.print(lli+ ", ") ;
        }
        System.out.println();
        
        System.out.println("KIJ matrix initialized") ;
       
    }
    
}
