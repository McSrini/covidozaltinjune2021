/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ca.mcmaster.covidozaltinjune2021.drivers;

import static ca.mcmaster.covidozaltinjune2021.Constants.*;
import ca.mcmaster.covidozaltinjune2021.cplex.*;
 

/**
 *
 * @author tamvadss
 */
public class Driver {
    
    public static void main(String[] args) throws Exception{
         
        Base_ModelBuilder builder_up = new Upper_ModelBuilder ();
        builder_up.solve( );
        
        Base_ModelBuilder builder_low = new Lower_ModelBuilder ();
        builder_low.solve( );
         
    }
    
}
