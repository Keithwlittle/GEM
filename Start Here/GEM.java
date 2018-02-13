import java.util.StringTokenizer;
import java.io.*;
import java.util.Random;//KWL 3/8/12 change
//KWL 3/22/12 change -- global replace MacCormack with MacCormack
public class GEM {
    static String dummystring,externalAfile,externalBfile,externalRVfile;
    static boolean EquationSolver, Linear, Dynamic,DynParm,Rnonlinear,AnalDeriv,SaveNewton,NegProblem,OutputTime,Explicit;
    static boolean CallFlowandEMultipliers, CallWiggle;
    static int NEQN, MaxNSV,MaxNSVxNT,NI, NB, NT,MAS,MixLengthOption, LGBCOption, WiggleOption,MaxNewtonIteration;
    static int MaxSamplingIteration, NumSamples;
    static int NColCompartment, NColPerInface, NColPerEcoef, NColPerFlow, NColKdandTempCoef, NColBoundary, NColInitial;
    static int DynType, NumTimeSteps,OutputTimeStep, DynPrint,NumTimePrint,NumSectnPrint,NumPrintInterval,VerboseLevel;
    static float Alpha, FlowTolerance,DelT,DeltFD,UB,LB,SC;
    static int SVComp[][], NEQNsm[], EqnToSectn[][], SectnToEqn[][];
    static double A[][], b[], c[], BBC[], V[], R[], RVEqnSolver[], cupdate[], f[],G[],XJAC[][];
    static float Compartment[][], Inface[][], ECoef[][], KdandTempCoef[][], Boundary[][], Flow[][], AlphaWeighted[][];
    static double cinitial[][],cinitialext[];
    static double AC[], RVt[], RVtp1[];
    static float Asm[][], Bsm[], MixLength;
    static double Aimplicit[][],Bimplicit[];
    static double Slope1[],Slope2[],chat[];
    //static double Asolve[][]; KWL 3/5/12 change
    static double DTmax; // max time step for stability
    static double StableCriterion;
    static double NewtonSF=1.0;
    static double Epsilon;
    static double xTrial[],NewtStep[];
    static double RVC[],RVCtp1[],ACFB[],ACFBtp1[];//KWL 3/9/12 change
    static boolean Shellfofc,Shelldfdc, ShellRofc, ShelldRdc,ShellCompartment,ShellInterfaces,ShellFlows,ShellECoefficients,ShellFlowandEMultipliers;
    static boolean ShellLinearKdandTempCoef,ShellBoundary,ShellLoads,ShellAEqnSolver,ShellBEqnSolver,ShellRVEqnSolver;
    //MB code
    static boolean DoMassBalance;
    static double MBBoundary[][];
    static double MBExtVolSinks[][];
    static double MBExtAreaSinks[][];
    static double SumMassInSV[];
    static double SumMassOutSV[][];
    static double MBInitial[];
    static double MBResidual[];
    static double MBLoad[];
    static int MBBoundaryLength; //length of MBBoundary array
    static int MBBoundaryRowCount;// will count rows of MBBoundary array
    static int MBExtVolSinksRowCount;// will count rows of MBExtVolSinks array
    static int MBExtAreaSinksRowCount;//will count rows of MBExtAreaSinks array
 
    static double BoundaryTransportFlux[][];
    static double BoundaryArealFlux[][];
    
  
      public static void main(String args[]) {
      	      VerboseLevel = 1;//1 is minimal System.out 2 is next, 3 is most
       	      if(VerboseLevel > 1) System.out.println("I'm in class main() method");
		readControlParms();
		if(VerboseLevel == 3)System.out.println("EquationSolver = " + EquationSolver);
		if(EquationSolver == false) readSVCompmap();
		ArrayDeclaration();

                // Start solutions
		if (Linear == false && Explicit == false ) {
			NonlinearInitialize();
		}

		if (Dynamic == false) {
			SteadyStateOption();
			System.out.println("  We're done with the steady-state simulation!");
 		}
		else { // dynamic solution
			GetReady();
			TimeLoop();
			System.out.println("  We're done with the dynamic simulation!");
		}// end if Dynamic == false

	 }// end Main

	static void readControlParms(){
		try
		{
			if(VerboseLevel > 1) System.out.println("I'm in readControlParms method");
			String strFile = "Control.csv";
			BufferedReader br = new BufferedReader( new FileReader(strFile));//create BufferedReader to read csv file
			String strLine = "";
			StringTokenizer st = null;
			strLine = br.readLine();//'Environmental System or Equation Solver mode (EnvS, EqnS)?
			if(VerboseLevel == 3)	System.out.println (strLine);
			dummystring = br.readLine();//EnvS,EqnS
			dummystring = dummystring.toUpperCase();
			if(VerboseLevel == 3)	System.out.println (dummystring);
			if (dummystring.equals("EQNS")) {
				EquationSolver = true;
			}
			else {
				if (dummystring.equals("ENVS")) {
					EquationSolver = false;
				}
				else {
					System.out.println("  Error: neither EQNS or ENVS");
					System.exit(1);
				}
			}
			if(VerboseLevel == 3)	System.out.println ("EquationSolver is " + EquationSolver);

			strLine = br.readLine();//''Linear problem? (Y,N)
			if(VerboseLevel == 3)	System.out.println (strLine);
	
			dummystring = br.readLine();
			dummystring = dummystring.toUpperCase();
			if(VerboseLevel == 3)	System.out.println (dummystring);

			if (dummystring.equals("Y")) {
				Linear = true;
			}
			else {
				if (dummystring.equals("N")) {
					Linear = false;
				}
				else {
					System.out.println("  Error: Linear neither Y or N");
					System.exit(1);
				}
			}
			if(VerboseLevel == 3)	System.out.println ("Linear is " + Linear);

			strLine = br.readLine();//''Steady-state simulation (else dynamic)? (Y,N)
			if(VerboseLevel == 3)	System.out.println (strLine);
	
			dummystring = br.readLine();
			dummystring = dummystring.toUpperCase();
			if(VerboseLevel == 3)	System.out.println (dummystring);
			
			if (dummystring.equals("Y")) {
				Dynamic = false;
			}
			else {
				if (dummystring.equals("N")) {
					Dynamic = true;
				}
				else {
					System.out.println("  Error: Dynamic neither Y nor N");
					System.exit(1);
				}
			}
			if(VerboseLevel == 3)	System.out.println ("Dynamic is " + Dynamic);
		

			strLine = br.readLine();//''If EqnS, enter the number of equations.
			if(VerboseLevel == 3){
				System.out.println (strLine);
			}			
			dummystring = br.readLine();
			NEQN = Integer.parseInt(dummystring);
			if(VerboseLevel == 3){
				System.out.println ("NEQN is " + NEQN);
			}			

			strLine = br.readLine();//''Maximum number of state variables, number of nonboundary compartments, number of dummy + boundary compartments
			if(VerboseLevel == 3){
				System.out.println (strLine);
			}			
			strLine = br.readLine();//read comma-separated MaxNSV,NI,NB text line
			st = new StringTokenizer(strLine, ",");
			MaxNSV = Integer.parseInt(st.nextToken());
			if(VerboseLevel == 3){
				System.out.println ("MaxNSV = "+ MaxNSV);
			}
			NI = Integer.parseInt(st.nextToken());
			if(VerboseLevel == 3){
				System.out.println ("NI = "+ NI);
			}
			NB = Integer.parseInt(st.nextToken());
			if(VerboseLevel == 3){
				System.out.println ("NB = "+ NB);
			}

			strLine = br.readLine();//''Maximum number of adjacent compartments
			if(VerboseLevel == 3){
				System.out.println (strLine);
			}
			dummystring = br.readLine();
			MAS = Integer.parseInt(dummystring);
			if(VerboseLevel == 3){
				System.out.println ("MAS is " + MAS);
			}
			strLine = br.readLine();//If EnvS enter MixLengthOption
			if(VerboseLevel == 3){
				System.out.println (strLine);
			}
			dummystring = br.readLine();
			MixLengthOption = Integer.parseInt(dummystring);
			if(VerboseLevel == 3){
				System.out.println ("MixLengthOption = "+ MixLengthOption);
			}

			strLine = br.readLine();//''Enter the around-compartments flow balance tolerance (%).
			if(VerboseLevel == 3){
				System.out.println (strLine);
			}
			dummystring = br.readLine();
			FlowTolerance = Float.parseFloat(dummystring);
			if(VerboseLevel == 3){
				System.out.println ("FlowTolerance is " + FlowTolerance);
			}
			strLine = br.readLine();//''If dynamic, enter 1 for FT,Euler's method, 2 for CT,MacCormack's method, 3 for BT,Implicit method, 4 for CT,Crank-Nicholson
			if(VerboseLevel == 3){
				System.out.println (strLine);
			}
			dummystring = br.readLine();
			DynType = Integer.parseInt(dummystring);
			if(VerboseLevel == 3){
				System.out.println ("DynType is " + DynType);
			}
			if (DynType > 4) {
				System.out.println("  Error: DynType > 4");
				System.exit(1);
			}
			strLine = br.readLine();//'If dynamic, enter the time step in days.
			if(VerboseLevel == 3){
				System.out.println (strLine);
			}
			dummystring = br.readLine();
			DelT = Float.parseFloat(dummystring);
			if(VerboseLevel == 3){
				System.out.println ("DelT is " + DelT);
			}
			strLine = br.readLine();//''If dynamic, enter the number of time steps.
			if(VerboseLevel == 3){
				System.out.println (strLine);
			}
			dummystring = br.readLine();
			NumTimeSteps = Integer.parseInt(dummystring);
			if(VerboseLevel == 3){
				System.out.println ("NumTimeSteps is " + NumTimeSteps);
			}
			strLine = br.readLine();//''if dynamic, enter the time step at which to output diagnostic files
			if(VerboseLevel == 3){
				System.out.println (strLine);
			}
			dummystring = br.readLine();
			OutputTimeStep = Integer.parseInt(dummystring);
			if(VerboseLevel == 3){
				System.out.println ("OutputTimeStep is " + OutputTimeStep);
			}
			strLine = br.readLine();//'If dynamic and EnvS, enter option (1,2, or 3) for writing outputs.
			if(VerboseLevel == 3){
				System.out.println (strLine);
			}
			dummystring = br.readLine();
			DynPrint = Integer.parseInt(dummystring);
			if(VerboseLevel == 3){
				System.out.println ("DynPrint is " + DynPrint);
			}
			if (DynPrint > 4) {
				System.out.println("  Error: DynPrint > 4");
				System.exit(1);
			}

			strLine = br.readLine();//''If write option = 1, enter the time step at which to print the profile.
			if(VerboseLevel == 3){
				System.out.println (strLine);
			}
			dummystring = br.readLine();
			NumTimePrint = Integer.parseInt(dummystring);
			if(VerboseLevel == 3){
				System.out.println ("NumTimePrint is " + NumTimePrint);
			}
			strLine = br.readLine();//if write option = 2, enter the compartment number and the time step print interval.
			if(VerboseLevel == 3){
				System.out.println (strLine);
			}
			strLine = br.readLine();
			st = new StringTokenizer(strLine, ",");
			NumSectnPrint = Integer.parseInt(st.nextToken());
			if(VerboseLevel == 3){
				System.out.println ("NumSectnPrint = "+ NumSectnPrint);
			}
			NumPrintInterval = Integer.parseInt(st.nextToken());
			if(VerboseLevel == 3){
				System.out.println ("NumPrintInterval = "+ NumPrintInterval);
			}
			strLine = br.readLine();//''If nonlinear and dynamic, is your R matrix is nonlinear, i.e. a function of the concentrations? (Y,N, N/A)
			if(VerboseLevel == 3){
				System.out.println (strLine);
			}
			dummystring = br.readLine();
			dummystring = dummystring.toUpperCase();
			if(VerboseLevel == 3){
				System.out.println (dummystring);
			}
			Rnonlinear = false;//default
			if (dummystring.equals("Y")) {
				Rnonlinear = true;
			}
			if(VerboseLevel == 3){
				System.out.println ("Rnonlinear is " + Rnonlinear);
			}
			strLine = br.readLine();//''If nonlinear, are analytical derivatives for the Jacobian available (Y,N, N/A)?
			if(VerboseLevel == 3){
				System.out.println (strLine);
			}
			dummystring = br.readLine();
			dummystring = dummystring.toUpperCase();
			if(VerboseLevel == 3){
				System.out.println (dummystring);
			}
			AnalDeriv = false;// default
			if (dummystring.equals("Y")) {
				AnalDeriv = true;
			}
			if(VerboseLevel == 3){
				System.out.println ("AnalDeriv is " + AnalDeriv);
			}
			strLine = br.readLine();//'If no analytical derivatives, enter the finite difference multiplier for derivative estimation
			if(VerboseLevel == 3){
				System.out.println (strLine);
			}
			dummystring = br.readLine();
			DeltFD = Float.parseFloat(dummystring);
			if(VerboseLevel == 3){
				System.out.println ("DeltFD is " + DeltFD);
			}
			strLine = br.readLine();//''If nonlinear, would you like to save Newton Algorithm diagnostics to NewtonInfo.txt file?
			if(VerboseLevel == 3){
				System.out.println (strLine);
			}
			dummystring = br.readLine();
			dummystring = dummystring.toUpperCase();
			if(VerboseLevel == 3){
				System.out.println (dummystring);
			}
			SaveNewton = false; //default
			if (dummystring.equals("Y")) {
				SaveNewton = true;
			}
			if(VerboseLevel == 3){
				System.out.println ("SaveNewton is " + SaveNewton);
			}
			strLine = br.readLine();//''If nonlinear, will negative concentrations present numerical problems during Newton method iterations?
			if(VerboseLevel == 3){
				System.out.println (strLine);
			}
			dummystring = br.readLine();
			dummystring = dummystring.toUpperCase();
			if(VerboseLevel == 3){
				System.out.println (dummystring);
			}
			NegProblem = false; //default
			if (dummystring.equals("Y")) {
				NegProblem = true;
			}
			if(VerboseLevel == 3){
				System.out.println ("NegProblem is " + NegProblem);
			}
			strLine = br.readLine();//''If nonlinear, enter a reasonable minimum and maximum value for your concentrations.
			if(VerboseLevel == 3){
				System.out.println (strLine);
			}
			strLine = br.readLine();
			st = new StringTokenizer(strLine, ",");
			LB = Float.parseFloat(st.nextToken());
			if(VerboseLevel == 3){
				System.out.println ("LB = "+ LB);
			}
			UB = Float.parseFloat(st.nextToken());
			if(VerboseLevel == 3){
				System.out.println ("UB = "+ UB);
			}
			strLine = br.readLine();//''If nonlinear, Newton's Method stopping criterion
			if(VerboseLevel == 3){
				System.out.println (strLine);
			}
			dummystring = br.readLine();
			SC = Float.parseFloat(dummystring);
			if(VerboseLevel == 3){
				System.out.println ("SC is " + SC);
			}
			strLine = br.readLine();//iIf nonlinear, enter (1) MaxNewtonIteration, (2) MaxSamplingIteration, (3) NumSamples
			if(VerboseLevel == 3){
				System.out.println (strLine);
			}
			strLine = br.readLine();
			st = new StringTokenizer(strLine, ",");
			MaxNewtonIteration = Integer.parseInt(st.nextToken());
			if(VerboseLevel == 3){
				System.out.println ("MaxNewtonIteration = "+ MaxNewtonIteration);
			}
			MaxSamplingIteration = Integer.parseInt(st.nextToken());
			if(VerboseLevel == 3){
				System.out.println ("MaxSamplingIteration = "+ MaxSamplingIteration);
			}
			NumSamples = Integer.parseInt(st.nextToken());
			if(VerboseLevel == 3){
				System.out.println ("NumSamples = "+ NumSamples);
			}
			strLine = br.readLine();//If nonlinear, are you providing nonlinear f(c) source/sink inputs with user-provided fofc.exe and the Shell functionality? (Y/N)
			if(VerboseLevel == 3){
				System.out.println (strLine);
			}
			dummystring = br.readLine();
			dummystring = dummystring.toUpperCase();
			if(VerboseLevel == 3){
				System.out.println (dummystring);
			}
			Shelldfdc = false; //default
			if (dummystring.equals("Y")) {
				Shellfofc = true;
				if(AnalDeriv == true)Shelldfdc = true;//you've got to provide the derivatives via Shell
			}else{
				Shellfofc = false;
			}
			if(VerboseLevel == 3) System.out.println("Shellfofc = "+Shellfofc);
			strLine = br.readLine();//if nonlinear and dynamic, are you providing nonlinear R(c)inputs with user-provided Rofc.exe and the Shell functionality? (Y/N)
			if(VerboseLevel == 3){
				System.out.println (strLine);
			}
			dummystring = br.readLine();
			dummystring = dummystring.toUpperCase();
			if(VerboseLevel == 3){
				System.out.println (dummystring);
			}
			ShelldRdc = false; //default
			if (dummystring.equals("Y")) {
				ShellRofc = true;
				if(AnalDeriv == true)ShelldRdc = true;//you've got to provide the derivatives via Shell
			}else{
				ShellRofc = false;
			}
			if(VerboseLevel == 3) System.out.println("ShellRofc = "+ShellRofc);
			
			strLine = br.readLine();//If dynamic and EnvS, do you want to use the Shell functionality to update Compartment.csv? (Y,N, N/A)
			if(VerboseLevel == 3){
				System.out.println (strLine);
			}
			dummystring = br.readLine();
			dummystring = dummystring.toUpperCase();
			if(VerboseLevel == 3){
				System.out.println (dummystring);
			}
			if (dummystring.equals("Y")) {
				ShellCompartment = true;
				DynParm = true;
			}else{
				ShellCompartment = false;
			}
			if(VerboseLevel == 3) System.out.println("ShellCompartment = "+ShellCompartment);
			
			strLine = br.readLine();//If dynamic and EnvS, do you want to use the Shell functionality to update Interfaces.csv? (Y,N, N/A)
			if(VerboseLevel == 3){
				System.out.println (strLine);
			}
			dummystring = br.readLine();
			dummystring = dummystring.toUpperCase();
			if(VerboseLevel == 3){
				System.out.println (dummystring);
			}
			if (dummystring.equals("Y")) {
				ShellInterfaces = true;
				DynParm = true;
			}else{
				ShellInterfaces = false;
			}
			if(VerboseLevel == 3) System.out.println("ShellInterfaces = "+ShellInterfaces);
			
			strLine = br.readLine();//If dynamic and EnvS, do you want to use the Shell functionality to update Flows.csv? (Y,N, N/A)
			if(VerboseLevel == 3){
				System.out.println (strLine);
			}
			dummystring = br.readLine();
			dummystring = dummystring.toUpperCase();
			if(VerboseLevel == 3){
				System.out.println (dummystring);
			}
			if (dummystring.equals("Y")) {
				ShellFlows = true;
				DynParm = true;
			}else{
				ShellFlows = false;
			}
			if(VerboseLevel == 3) System.out.println("ShellFlows = "+ShellFlows);
 
			strLine = br.readLine();//If dynamic and EnvS, do you want to use the Shell functionality to update ECoefficients.csv? (Y,N, N/A)
			if(VerboseLevel == 3){
				System.out.println (strLine);
			}
			dummystring = br.readLine();
			dummystring = dummystring.toUpperCase();
			if(VerboseLevel == 3){
				System.out.println (dummystring);
			}
			if (dummystring.equals("Y")) {
				ShellECoefficients = true;
				DynParm = true;
			}else{
				ShellECoefficients = false;
			}
			if(VerboseLevel == 3) System.out.println("ShellECoefficients = "+ShellECoefficients);
 
			strLine = br.readLine();//If dynamic and EnvS, do you want to use the Shell functionality to update FlowandEMultipliers.csv? (Y,N, N/A)
			if(VerboseLevel == 3){
				System.out.println (strLine);
			}
			dummystring = br.readLine();
			dummystring = dummystring.toUpperCase();
			if(VerboseLevel == 3){
				System.out.println (dummystring);
			}
			if (dummystring.equals("Y")) {
				ShellFlowandEMultipliers = true;
				DynParm = true;
			}else{
				ShellFlowandEMultipliers = false;
			}
			if(VerboseLevel == 3) System.out.println("ShellFlowandEMultipliers = "+ShellFlowandEMultipliers);
			
			strLine = br.readLine();//If dynamic and EnvS, do you want to use the Shell functionality to update LinearKdandTempCoef.csv? (Y,N, N/A)
			if(VerboseLevel == 3){
				System.out.println (strLine);
			}
			dummystring = br.readLine();
			dummystring = dummystring.toUpperCase();
			if(VerboseLevel == 3){
				System.out.println (dummystring);
			}
			if (dummystring.equals("Y")) {
				ShellLinearKdandTempCoef = true;
				DynParm = true;
			}else{
				ShellLinearKdandTempCoef = false;
			}
			if(VerboseLevel == 3) System.out.println("ShellLinearKdandTempCoef = "+ShellLinearKdandTempCoef);
			
			strLine = br.readLine();//If dynamic and EnvS, do you want to use the Shell functionality to update Boundary.csv? (Y,N, N/A)
			if(VerboseLevel == 3){
				System.out.println (strLine);
			}
			dummystring = br.readLine();
			dummystring = dummystring.toUpperCase();
			if(VerboseLevel == 3){
				System.out.println (dummystring);
			}
			if (dummystring.equals("Y")) {
				ShellBoundary = true;
				DynParm = true;
			}else{
				ShellBoundary = false;
			}
			if(VerboseLevel == 3) System.out.println("ShellBoundary = "+ShellBoundary);
			
			strLine = br.readLine();//If dynamic and EnvS, do you want to use the Shell functionality to update Loads.csv? (Y,N, N/A)
			if(VerboseLevel == 3){
				System.out.println (strLine);
			}
			dummystring = br.readLine();
			dummystring = dummystring.toUpperCase();
			if(VerboseLevel == 3){
				System.out.println (dummystring);
			}
			if (dummystring.equals("Y")) {
				ShellLoads = true;
				DynParm = true;
			}else{
				ShellLoads = false;
			}
			if(VerboseLevel == 3) System.out.println("ShellLoads = "+ShellLoads);
			
			strLine = br.readLine();//If dynamic and EnvS, do you want to use the Shell functionality to update AEqnSolver.csv? (Y,N, N/A)
			if(VerboseLevel == 3){
				System.out.println (strLine);
			}
			dummystring = br.readLine();
			dummystring = dummystring.toUpperCase();
			if(VerboseLevel == 3){
				System.out.println (dummystring);
			}
			if (dummystring.equals("Y")) {
				ShellAEqnSolver = true;
				DynParm = true;
			}else{
				ShellAEqnSolver = false;
			}
			if(VerboseLevel == 3) System.out.println("ShellAEqnSolver = "+ShellAEqnSolver);
			
			strLine = br.readLine();//If dynamic and EnvS, do you want to use the Shell functionality to update BEqnSolver.csv? (Y,N, N/A)
			if(VerboseLevel == 3){
				System.out.println (strLine);
			}
			dummystring = br.readLine();
			dummystring = dummystring.toUpperCase();
			if(VerboseLevel == 3){
				System.out.println (dummystring);
			}
			if (dummystring.equals("Y")) {
				ShellBEqnSolver = true;
				DynParm = true;
			}else{
				ShellBEqnSolver = false;
			}
			if(VerboseLevel == 3) System.out.println("ShellBEqnSolver = "+ShellBEqnSolver);
			
			strLine = br.readLine();//If dynamic and EnvS, do you want to use the Shell functionality to update RVEqnSolver.csv? (Y,N, N/A)
			if(VerboseLevel == 3){
				System.out.println (strLine);
			}
			dummystring = br.readLine();
			dummystring = dummystring.toUpperCase();
			if(VerboseLevel == 3){
				System.out.println (dummystring);
			}
			if (dummystring.equals("Y")) {
				ShellRVEqnSolver = true;
				DynParm = true;
			}else{
				ShellRVEqnSolver = false;
			}
			if(VerboseLevel == 3) System.out.println("ShellRVEqnSolver = "+ShellRVEqnSolver);
	 
    		br.close(); // close input stream
		}// end try
		catch(Exception e)
		{
			System.out.println("Exception while reading Control.csv file: " + e);
		}// end catch

		//DynParm = false; // !! hard-wired 
	} //end readControlParms method

	static void readSVCompmap(){
		if(VerboseLevel > 1){
				System.out.println ("I'm in readSVCompmap method");
			}
		NT = NI + NB;
		MaxNSVxNT = MaxNSV*NT;
		NEQN = 0;
		if(VerboseLevel == 3){
				System.out.println ("MaxNSV = " + MaxNSV);
				System.out.println("MaxNSVxNT = " + MaxNSVxNT);		
			}

		SVComp = new int[MaxNSV+1][NT+1];// create new 2-D array object

		try{
			String strFile = "SVCompmap.csv";
			BufferedReader br = new BufferedReader( new FileReader(strFile));//create BufferedReader to read csv file
			String strLine = "";
			StringTokenizer st = null;
			strLine = br.readLine();//'"SV","Comp 1","Comp 2","…",,,,"Comp NT"
			if(VerboseLevel == 3){
				System.out.println (strLine);
			}
			int i;
			int j;
			int dummySV;
			for (i=1; i < MaxNSV+1;i++) {
				if(VerboseLevel == 3){
				System.out.println ("i = " + i);
			}
				strLine = br.readLine();
				if(VerboseLevel == 3){
				System.out.println (strLine);
			}
				st = new StringTokenizer(strLine, ",");
				dummySV = Integer.parseInt(st.nextToken());
				if(VerboseLevel == 3){
				System.out.println ("dummySV = " + dummySV);
			}
				for (j = 1; j < NT+1 ;j++) {
					if(VerboseLevel == 3){
						System.out.println ("j = " + j);
					}
					SVComp[i][j] = Integer.parseInt(st.nextToken());
					if(VerboseLevel == 3){
						System.out.println ("SVComp["+i+"]["+j+"] = " +SVComp[i][j]);
					}
					if (SVComp[i][j] > 1) {
						System.out.println("  Error: SVComp[][] element > 1");
						System.exit(1);
				}else {
					NEQN = NEQN+SVComp[i][j];
					if(VerboseLevel == 3){
						System.out.println ("NEQN = " + NEQN);
					}
					}
				} // end j loop
			} // end i loop
			br.close();
		} // end try
		catch(Exception e)
		{
			System.out.println("Exception while reading Control.csv file: " + e);
		} // end catch
	} // end readSVCompmap method


	static void ArrayDeclaration(){
		if(VerboseLevel > 1){
			System.out.println( "I'm in ArrayDeclaration()");
		}
		//EnvS and EqnS ...
		A = new double[NEQN+1][NEQN+1];// create new 2-D array object
		b = new double[NEQN+1];
		c = new double[NEQN+1];
		Explicit = false; // default
		if (Dynamic == true) {
			if (DynType ==1) {
				if (Rnonlinear == false) {
					Explicit = true;
				}
			}
		}
		if(VerboseLevel == 3){
			System.out.println("Explicit =  " + Explicit);
		}

		if (Linear == false) {
			f = new double[NEQN+1];
			if (Explicit == false) {//nonlinear and Newton's method
				
				if(Linear == false){
				f = new double[NEQN+1];
				if(Explicit == false){//nonlinear and Newton's method
					if(Dynamic == true && DynType == 2){
						System.out.println("MacCormack method inappropriate for nonlinear systems");
						System.exit(1);
					}
				}// end if(Explicit == false)
				G = new double[NEQN+1];
				XJAC = new double[NEQN+1][NEQN+1];
				//System.out.println("VB code to open NewtonInfo.dng not yet translated in ArrayDeclaration !!");//KWL 3/6/12 change
				}//end if (Linear == false)
			}// end if Explicit == false
		}// end if Linear == false

		if (EquationSolver == true) {
   			externalAfile = "AEqnSolver.csv";
   			if(VerboseLevel == 3) System.out.println("externalAfile = "+externalAfile);
   			externalBfile = "BEqnSolver.csv";
   			if(VerboseLevel == 3) System.out.println("externalBfile = "+externalBfile);
   			if(Dynamic == true){
   				if(DynType == 2){
   					System.out.println("Error: MacCormack Method Not Applicable to Equation Solver Mode");
   					System.exit(1);
   				}
   				RVEqnSolver = new double[NEQN+1]; //Diagonal elements in R*V matrix in (dRVc/dt) = Ac + b
   				String ExternalInitialFile = "InitialEqnSolver.csv";
   			}
 		}
		else { // EquationSolver = false
			if(VerboseLevel > 1){
				System.out.println("I'm in ArrayDeclaration");
			}
			NEQNsm = new int[MaxNSV+1];
			EqnToSectn = new int[MaxNSV+1][NT+1];
			SectnToEqn = new int[MaxNSV+1][NT+1];
			BBC = new double[NEQN+1];
			V = new double[NEQN+1];
			R = new double[NEQN+1];
    			NColCompartment = 10; // number of columns in array Compartment
    			NColPerInface = 5; //'number of columns in array Inface for each interface
    			NColPerEcoef = 2; //' number of columms in array E for each interface
    			NColPerFlow = 2; // number of columms in array Q for each interface
    			NColKdandTempCoef = 4; // number of columns in array KdandTempCoef
    			NColBoundary = 3; // number of columns in array Boundary
    			NColInitial = 3; //number of columms in cinitial array
			Compartment = new float[NT+1][NColCompartment+1];
			Inface = new float[NT+1][2+NColPerInface*MAS];
			ECoef = new float[MaxNSVxNT+1][3 + NColPerEcoef * MAS];
			KdandTempCoef = new float[MaxNSVxNT+1][NColKdandTempCoef+1];
			Boundary= new float[MaxNSVxNT+1][NColBoundary+1];
			Flow = new float[MaxNSVxNT+1][3 + NColPerFlow * MAS];
			AlphaWeighted = new float[NT+1][2*MAS + 2];
			// following line changed by KWL 3/5/12
			CallFlowandEMultipliers = true;//hard-wired by user to false to speed execution if sure no flow and e multipliers are applicable
			CallWiggle = false;//initialize
		}
	} // end ArrayDeclaration method

	static void NonlinearInitialize(){
		int II,JJ;
		double x;
	       int Count = 0; //counts equations (to NEQN)
		int Count1 = 0; // counts rows of MaxNSVxNT arrays
                cupdate= new double[NEQN+1];
                // determine machine epsilon
                Epsilon = 1;
                while (1+Epsilon>1) {
                    Epsilon = Epsilon/2;
                }
              if (VerboseLevel == 3) System.out.println("machine epsilon = "+Epsilon);
              if (EquationSolver == true) {
                	try{
                		String strFile = "GuessEqnSolver.csv";
                		BufferedReader br = new BufferedReader( new FileReader(strFile));//create BufferedReader to read csv file
                		String strLine = "";
                		for (int i=1; i < NEQN+1;i++) {
                			strLine = br.readLine();
                			c[i] = Double.parseDouble(strLine);//initial guess
                			c[i]=c[i]*NewtonSF;
                			if(VerboseLevel == 3) {
                				System.out.println("NewtonSF = "+NewtonSF);
                				System.out.println("c["+i+"] = "+c[i]);
                			}
                		}
                	} // end try
                	catch(Exception e)
                	{
                            System.out.println("Exception while reading GuessEqnSolver.csv file: " + e);
                        } // end catch
            } else{// EquationSolver false
                       try{
                    	if(VerboseLevel> 1){
					System.out.println("I'm in nonlinear initialize method, equationsolver false");
				}    
			String strFile = "Guess.csv";
			BufferedReader br = new BufferedReader( new FileReader(strFile));//create BufferedReader to read csv file
			String strLine = "";
			StringTokenizer st = null;
			strLine = br.readLine();//SV,Compartment,Initial Guess(g/m3)
			if(VerboseLevel == 3){
				System.out.println(strLine);
			}
			for (int i=1; i < MaxNSV+1;i++) {
                            for (int j=1; j<NT+1; j++) {
                                Count1++;
                                strLine = br.readLine();
                                if(VerboseLevel == 3){
                                	System.out.println(strLine);
                                }
				st = new StringTokenizer(strLine, ",");
				II = Integer.parseInt(st.nextToken());
				JJ = Integer.parseInt(st.nextToken());
				x = Double.parseDouble(st.nextToken());
				//System.out.println("II = "+II+" JJ = "+JJ+" x = "+x);//SV, Compartment, guess for c
				if(Compartment[j][3] == 0){
					if(SVComp[i][j] == 1){
						Count++;
						if(Dynamic == false){
							c[Count]=x*NewtonSF;
							if(c[Count] == -999) {
									System.out.println("Error: initial guess is -999 in NonLinearInitialize");
									System.exit(1);
							}
						}else{
							cupdate[Count] = x*NewtonSF;
							if(cupdate[Count] == -999) {
									System.out.println("Error: initial guess is -999 in NonLinearInitialize");
									System.exit(1);
							}
						}
					}//end if(SVComp[i][j] == 1)
				}//if(Compartment[j][3] == 0)
                            }
			} // end i loop
			br.close();
                    } // end try
                    catch(Exception e)
                    {
                            System.out.println("Exception while reading Guess.csv file: " + e);
                    } // end catch
               
			if(Count != NEQN){
				System.out.println("Error: Count != NEQN in NonLinearInitialize()");
				System.exit(1);
			}
			if(Count1 != MaxNSVxNT){
				System.out.println("Error: Count1 != MaxNSVxNT in NonLinearInitialize()");
				System.exit(1);
			}
            }//end if(EquationSolver == true)
	} // end NonlinearInitialize method


    static void SteadyStateOption() {
    	    int I,J;
    	    //Asolve = new double[NEQN+1][NEQN+2];//will hold RHS as last column for LinearSolver() method -- KWL 3/5/12 change
    	    ParmUpdate(0);
	
	if(EquationSolver == false){
		//load transport terms into array A
		LoadATrans();//BBC loaded with boundary conditions only and on LHS
		LoadALinearSrcSnk();
		if(CallWiggle == true) Wiggle(0);
		ForceFn(0);
	}else{
		ForceFn(0);
	}//end if(EquationSolver == false)
	Aout(0);// KWL 3/14 change
	Bout(0);// KWL 3/14 change
	if (Linear == true){
		//start KWL 3/5/12 change
		// load Asolve matrix for solve method putting b vector back to RHS
		//for(I=1;I<NEQN+1;I++){
			//for(J=1;J<NEQN+2;J++){
				//if(J<NEQN+1){
					//Asolve[I][J]=A[I][J];
				//}else{
					//Asolve[I][J]=-b[I];//negative to get back to RHS
				//}
			//}
		//}
		//if(VerboseLevel == 3){
				//System.out.println("before LinearSolver, Asolve(NEQN,NEQN)= "+Asolve[NEQN][NEQN]);
		//}
		//LinearSolver(NEQN,Asolve,c);
		LinearSolver(NEQN,A,b);
		for (I=1;I<NEQN+1;I++){
			c[I] = b[I];
		}
		//end KWL 3/5/12 changes		
		if(VerboseLevel == 3){
			for (I = 1; I < NEQN+1; I++){
				System.out.println("c["+I+"]= "+c[I]);
			}	
		}
	}else{
		//KWL 3/15 change
		Newton(c,c,0);
	}
	ProfileSnapShot(-999);
	//KWL 3/15 change -- delete below VB code...
    }//end SteadyStateOption()
    
    static void GetReady() {
    	    
    	    	if(VerboseLevel> 1){
			System.out.println("I'm in GetReady()");
		}
		if (!EquationSolver) {
                        // create the cinitial array
                        cinitial = new double[MaxNSVxNT+1][NColInitial+1];
                         //READ IN initial.csv
                         try {
                         	 
                         	 String strFile = "Initial.csv";
                         	 BufferedReader br = new BufferedReader( new FileReader(strFile));//create BufferedReader to read csv file
                         	 String strLine = "";
                         	 StringTokenizer st = null;
                         	 strLine = br.readLine();//"SV","Compartment","Initial Condition"
                         	 if(VerboseLevel == 3){
                         	 	 System.out.println(strLine);
                         	 }
                         	 for (int i=1; i < MaxNSVxNT+1;i++) {
                         	 	 if(VerboseLevel == 3){
                         	 	 	 System.out.println("i = " + i);
                         	 	 }
                         	 	 strLine = br.readLine();
                         	 	 st = new StringTokenizer(strLine, ",");
                         	 	 for (int j = 1; j < NColInitial+1; j++){
                         	 	 	 cinitial[i][j] = Double.parseDouble(st.nextToken());
                         	 	 	 if(VerboseLevel == 3){
                         	 	 	 	 System.out.println("j = " + j);
                         	 	 	 	 System.out.println("cinitial = " +cinitial[i][j]);
                         	 	 	 }
                         	 	 } // end for j
                         	 } // end for i
                         	 br.close(); // close input stream
			 }// end try
			 catch(Exception e)
			 {
			 	System.out.println("Exception while reading initial.csv file: " + e);
			 }// end catch
		ParmUpdate(0);// initial (or only) update of parameters
                    } else { // EquationSolver is true
                        cinitialext= new double[NEQN+1];
                    } // end if !EquationSolver
        } // end GetReady() method
        
         static void TimeLoop(){
        	int TimeStep, ITP, I,J;
        	float TempTime = 0;
        	
        	OutputTime = false;
        	if (VerboseLevel > 1){
        		System.out.println("I'm in TimeLoop() and NumTimeSteps = " +NumTimeSteps);
        	}
        	float CumTime = 0;
        	
        	//MB code
        	DoMassBalance = true;//initialize
        	if(Dynamic == false) DoMassBalance = false;
        	//if(Linear == false) DoMassBalance = false;
        	if(DoMassBalance == true){
        		MBExtVolSinks = new double[MaxNSVxNT+1][3+1];//can't predict the length of this array, so just make it "big"
        		MBExtAreaSinks = new double[MaxNSVxNT+1][4+1];//can't predict the length of this array, so just make it "big"
        		SumMassInSV=new double[MaxNSV+1];//will cumulate mass entering modeled compartments
        		SumMassOutSV = new double[MaxNSV+1][3+1];//cumulate mass leaving modeled compartments, SumMassOutSV(*,1)will be transport, (*,2) will be volume-based sinks, (*,3) will be area-based sinks
        		MBInitial = new double[MaxNSV+1];
        		MBResidual = new double[MaxNSV+1];
        		MBLoad = new double[NEQN+1];
        		for(I=1;I<MaxNSV+1;I++){
        			SumMassInSV[I]=0;//initialize
        			for(J=1;J<4;J++){
        				SumMassOutSV[I][J] = 0;//initialize
        			}
        		}
        		//KWL 3/15 change -- delete VB code below ...
          	}// end if(DoMassBalance == true)
 
        	ITP = 0;
        	for (TimeStep = 1; TimeStep < NumTimeSteps+1; TimeStep++){
        		System.out.println ("  TimeStep = " +TimeStep);
        		/*
        		'Note -- TimeStep is equivalent to time step 't' in the documentation
        		'the update to C before the end of this loop updates C to 't+1'
        		'Therefore, any data loaded at time TimeStep will be 't'
        		'To load at 't+1', e.g. for BTCS RHS, need to pass TimeStep + 1 to loading subroutines
        		*/
        		ITP ++; // ITP is output print time counter
        		if (ITP > NumPrintInterval) {
        			ITP = 1;
        		}
        		if (TimeStep == OutputTimeStep) {
        			OutputTime = true;
        		}else{//KWL 3/15 change
        			OutputTime = false;
        		}
        		if (TimeStep == 1){ // initialize c vector, etc.
        			if (EquationSolver == true) {
        				for (I = 1; I < NEQN+1; I++){
        					c[I] = cinitialext[I];
        				}
        				TempTime = CumTime;
        			} else { // EquationSolver = false
        				int Count = 0;
        				int Count1 = 0;
        				for (I = 1; I < MaxNSV+1; I++){
        					for (J = 1; J < NT+1; J++){
        						Count1 ++;
        						if (Compartment[J][3] == 0 && SVComp[I][J] == 1) {
        							Count ++;
        							c[Count] = cinitial[Count1][NColInitial];
        							if (VerboseLevel==3){
        								System.out.println(" initial condition c[] ="  +c[Count]);
        							}
        							if (c[Count] == -999) {
        								System.out.println("Error: initial condition = -999 in TimeLoop method");
        								System.exit(1);
        							}
        						}
        					} // end J loop
        				} // end I loop
        				if (Count != NEQN) {
        					System.out.println("Error: Count != NEQN in TimeLoop method");
        					System.exit(1);
        				}
					//MB code
					//calculate MBBoundaryLength for mass balance 
					TransportFluxToBoundaryCompartments(0,0,0,0);//erase old and set up headers
					ArealFluxToBoundaryCompartments(0,0,0,0);//erase old and set up headers
					if(DoMassBalance == true){
						for(int ISV = 1;ISV < MaxNSV+1;ISV++){
							for(I = 1;I<NT+1;I++){
								if(Compartment[I][3] == 0){// interior compartment
									if(SVComp[ISV][I] == 1){
										int NAS = (int) Compartment[I][2];
										int ICInface = 2;
										for(J=1;J<NAS+1;J++){
											if(J == 1){
												ICInface = 2;	
											}else{
												ICInface += NColPerInface; 	
											}
											int IAS = (int) Inface[I][ICInface];
											if(Compartment[IAS][3] != 0){//boundary or dummy compartment
												MBBoundaryLength += 1;
											}
										}
									}
								}//end if(Compartment[I][3] == 0
							}// end for(I = 1;I<NT+1;I++){
						}//end for(int ISV = 1;ISV < MaxNSV+1;ISV++){
						MBBoundary = new double[MBBoundaryLength+1][5+1];// 5+1 silliness for consistency with VB code
					}// end if(DoMassBalance == true 					
            			} // end else EquationSolver = false
       			
        		} //end if TimeStep == 1

               // Select the cases
                    if (DynType == 1) {
                    	    EulerOption(TimeStep, CumTime,OutputTime,TempTime);
                    	    TempTime = CumTime + DelT;
                    }
                    if(DynType == 2){
                    	    MacCormack(TimeStep, CumTime,OutputTime,TempTime);
                    	    TempTime = CumTime + DelT;
                    }
 			if (DynType == 3) {
                            BackTimeOption(TimeStep,CumTime,OutputTime,TempTime);
                             TempTime = CumTime + DelT;
                    }

                    if (DynType == 4) {
                            CenterTimeOption(TimeStep,CumTime,OutputTime,TempTime);
                            TempTime = CumTime + DelT;
                    }
                    
                    if (GEM.DynPrint == 1){
                    	    if (TimeStep == NumTimePrint) {
                    	    	    ProfileSnapShot(TempTime);
                    	    }
                    }
                     if (GEM.DynPrint == 2){
                    	    if (ITP == 1) {
                    	    	    OneCompTimeSeries(NumSectnPrint,TempTime);
                    	    }
                    }                   
                    if (GEM.DynPrint == 3){
                    	    AllTimeSeries(TempTime);
                    }
                    
                    // MB code
                    if(DoMassBalance == true){
                    	    //update mass balance cumulators
                    	    //boundary transport first
                    	    int EqnBlkCounter,Eqn,Comp,AdjComp,ISV;
                    	    double CiCoef;//coefficient of "i" compartment
                    	    double CjCoef;//coefficient of "j" compartment
                    	    double BoundaryConc,Flux;
                    	    if (VerboseLevel == 3) System.out.println("in TimeLoop, MBBoundaryLength = "+MBBoundaryLength+" CumTime = "+CumTime);
                    	    for (I=1;I< MBBoundaryLength+1;I++){
                    	    	    ISV = (int) MBBoundary[I][1];
                    	    	    Comp = (int) MBBoundary[I][2];
                    	    	    AdjComp = (int) MBBoundary[I][3];
                    	    	    if (VerboseLevel == 3) System.out.println("in TimeLoop ISV,Comp,AdjComp = "+ISV+","+Comp+","+AdjComp);
                    	    	    CjCoef = MBBoundary[I][4];
                    	    	    CiCoef = MBBoundary[I][5];
				    if (VerboseLevel == 3) System.out.println("CjCoef = "+CjCoef+" CiCoef = "+CiCoef);
                    	    	    EqnBlkCounter = 0;
                    	    	    for(J=1;J<ISV-1+1;J++){
                    	    	    	  EqnBlkCounter += NEQNsm[J];  
                    	    	    }
                    	    	    Eqn = EqnBlkCounter+SectnToEqn[ISV][Comp];
                    	    	    // find boundary concentration
                    	    	    for(int K=1;K< MaxNSVxNT+1;K++){
                    	    	    	    if(Boundary[K][1] == ISV){//right SV
						    if(Boundary[K][2] == AdjComp) {// right boundary compartment
						    	    BoundaryConc = -999;//initialize
                                                            if(Compartment[AdjComp][3] == 1) {//dummy
                                                                    BoundaryConc = 0;// zeroes out CjCoef below
                                                            }
                                                            if(Compartment[AdjComp][3]==2){//fixed
                                                                   BoundaryConc = Boundary[K][3]; 
                                                            }
                                                            if(Compartment[AdjComp][3] == 3){//zero gradient
                                                                    BoundaryConc = 0;// 'redundant -- CjCoef already 0
                                                            }
                                                            if(Compartment[AdjComp][3] == 4){//linear gradient
                                                                    BoundaryConc = 0;//initialize
                                                                    int NAS = (int) Compartment[Comp][2];//NUMBER OF adjacent CompartmentS (INCLUDING BOUNDARY)
                                                                    //Get length (Avgl) between Comp and AdjComp
                                                                    float Avgl = GetLength(1, Comp, AdjComp);// use X,Y,Z coordinates
                                                                    //FIND NUMBER OF INTERIOR ADJ CompartmentS (NASINT) AND AVG LENGTH (AVGADJL)and concentration (cANB) OVER THESE
                                                                    int NASINT = 0;
                                                                    float AvgAdjL = 0;
                                                                    for(int IK=1;IK<NAS+1;IK++){
                                                                            int JC = 2;
                                                                            if(IK == 1){
                                                                                    JC = 2;
                                                                            }else{
                                                                                    JC += NColPerInface;
                                                                            }
                                                                            int IADJ = (int) Inface[Comp][JC];// ADJ Compartment NUMBER
                                                                            if (IADJ == -999){
                                                                                    System.exit(1);
                                                                                    System.out.println("Error: IADJ == -999 in TimeLoop");
                                                                            }
                                                                            if (Compartment[IADJ][3]== 0){//INTERIOR Compartment
                                                                                    NASINT+=1;
                                                                                    float AvglNew = GetLength(1,Comp,IADJ); //only X,Y,Z coordinate method used
                                                                                    AvgAdjL += AvglNew;
                                                                                    BoundaryConc += c[EqnBlkCounter+SectnToEqn[ISV][IADJ]];
                                                                            }
                                                                            if(NASINT > NAS){
                                                                                    System.exit(1);
                                                                                    System.out.println("Error: NASINT > NAS in TimeLoop");
                                                                            }
                                                                            if(NASINT ==0 ){
                                                                                    System.exit(1);
                                                                                    System.out.println("Error: (NASINT ==0  in TimeLoop");
                                                                            }
                                                                    }//end IK loop
                                                                    AvgAdjL = AvgAdjL/NASINT;
								    BoundaryConc = BoundaryConc/NASINT;//avg conc of adj interior compartments
								    BoundaryConc = -(Avgl / AvgAdjL) * BoundaryConc + (1 + Avgl / AvgAdjL) * c[Eqn];
                                                            }// end  if(Compartment[AdjComp][3] == 4){//linear gradient

                                                            if(BoundaryConc == -999){
                                                                    System.exit(1);
                                                                    System.out.println("Error: BoundaryConc == -999 in TimeLoop");
                                                            }
                                                        
                                                            Flux = (CiCoef*c[Eqn]+CjCoef*BoundaryConc)*DelT;//c is already Cd which is what is transported out
                                                            if (VerboseLevel == 3){
                                                                    System.out.println("in TimeLoop, CjCoef = "+CjCoef);
                                                                    System.out.println("in TimeLoop, BoundaryConc = "+BoundaryConc);
                                                                    System.out.println("in TimeLoop, CiCoef = "+CiCoef);
                                                                    System.out.println("in TimeLoop, c["+Eqn+"] = "+c[Eqn]);
                                                                    System.out.println("in TimeLoop,Flux = "+Flux);
                                                            }
                                                            if(Flux < 0){// net flux is from i to j, i.e. from compartment Comp to boundary compartment AdjComp
								    TempTime = CumTime+DelT;
                                                                    TransportFluxToBoundaryCompartments (TempTime,ISV,AdjComp,-Flux);
                                                                    SumMassOutSV[ISV][1] += -Flux;//transport
                                                            }else{
                                                                   SumMassInSV[ISV] += Flux; 
                                                            }
						    }// end if(Boundary[K][2] == AdjComp)
                    	    	    	    }//end  if(Boundary[K][1] == ISV)
                    	    	    }// end for(int K=1;K< MaxNSVxNT+1;K++)
                    	    }//end for (I=1;I< MBBoundaryLength+1;I++)
                    	    
                    	    //external loadings
                    	    Eqn = 0;
                    	    for(ISV=1;ISV<MaxNSV+1;ISV++){
                    	    	    for(I=1;I<NT+1;I++){
                    	    	    	    if(Compartment[I][3] == 0 && SVComp[ISV][I] == 1){//interior compartments only and SV ISV relevant to compartment I
                    	    	    	    	    Eqn += 1;
													//KWL change 9/26/13 to include dissolved and sorbed input load
                    	    	    	    	    //SumMassInSV[ISV] +=MBLoad[Eqn]*DelT/R[Eqn]; -- old way with dissolved only 
													SumMassInSV[ISV] +=MBLoad[Eqn]*DelT;
                    	    	    	    }
                    	    	    }
                    	    }
                    	    if(Eqn != NEQN){
                    	    	    System.out.println("Error: Eqn != NEQN in mass balance code for MBLoad");
                    	    	    System.exit(1);
                    	    }

    			// external volume-based source/sinks
    			for (I = 1;I<MaxNSVxNT+1;I++){
    				if(MBExtVolSinks[I][1] == -999){
    					break;// end of source/sink terms	
    				}else{
    					ISV = (int) MBExtVolSinks[I][1];
    					Comp = (int) MBExtVolSinks[I][2];
    					//first find equation corresponding to SV (MBExtVolSinks[1][*]) and compartment (MBExtVolSinks[2][*]
    					EqnBlkCounter = 0;
    					for(J=1;J<ISV-1+1;J++){//-1+1 silliness to agree with VB code
    						EqnBlkCounter+= NEQNsm[J];
    					}
    					Eqn = EqnBlkCounter+SectnToEqn[ISV][Comp];
    					if(MBExtVolSinks[I][3] < 0){// external sink
    						 SumMassOutSV[ISV][2]+=  Math.abs(MBExtVolSinks[I][3]) * c[Eqn] * DelT;
    					}else{//external source
    						SumMassInSV[ISV]+=  MBExtVolSinks[I][3] * c[Eqn] * DelT;
    					}
    				}
    			}//end (I = 1;I<MaxNSVxNT+1;I++)
    			//external area-based sinks (sources not currently allowed)
    			for(I=1;I<MaxNSVxNT+1;I++){
    				if(MBExtAreaSinks[I][1] == -999){
    					break;//end of source/sink terms
    				}else{
    					ISV = (int) MBExtAreaSinks[I][1];
    					int Comp1 = (int)MBExtAreaSinks[I][2];
    					int Comp2 = (int)MBExtAreaSinks[I][3];
					//first find equation corresponding to SV (MBExtAreaSinks(1,*)) and compartment (MBExtAreaSinks(2,*))
					EqnBlkCounter = 0;
					for(J=1;J<ISV-1+1;J++){//-1+1 silliness to agree with VB code
						EqnBlkCounter+= NEQNsm[J];
					}
					Eqn = EqnBlkCounter+SectnToEqn[ISV][Comp1];
					if(MBExtAreaSinks[I][4]<=0){// external sink
						Flux = Math.abs(MBExtAreaSinks[I][4])*c[Eqn]*DelT;
						TempTime = CumTime+DelT;
						ArealFluxToBoundaryCompartments (TempTime,ISV,Comp2,Flux);
						SumMassOutSV[ISV][3] +=Flux;
					}else{//external source
						System.out.println("Error in Mass Balance: external area-based source not allowed, user should specify as external load");
						System.exit(1);
					}
				}
    			}//end for(I=1;I<MaxNSVxNT+1;I++)
                   }//end if(DoMassBalance == true)                    
                     
                    CumTime +=DelT;
                    if (VerboseLevel== 3){
        		System.out.println("in TimeLoop(), CumTime = "+CumTime+" and TempTime = "+TempTime);
        	    }
                } // end TimeStep loop
                //MB code
                if(DoMassBalance == true){
                	
		
			//calculate initial mass in interior compartments
			int Count = 0;
			int Count1 = 0;
			for(int ISV=1;ISV<MaxNSV+1;ISV++){
				for(J=1;J<NT+1;J++){
					Count1 += 1;//will count rows of MaxNSVxNT arrays
					if(Compartment[J][3] == 0 && SVComp[ISV][J] == 1){
						Count += 1;
						MBInitial[ISV] += cinitial[Count1][NColInitial]*Compartment[J][7]*Compartment[J][9]*R[Count];
					}
				}
			}
			if(Count != NEQN){
				System.out.println("Error in TimeLoop: Count != NEQN");
				System.exit(1);
			}
			
                       //calculate residual mass in interior compartments
					   Count =0;//will count rows of MaxNSVxNT arrays
                       for(int ISV=1;ISV<MaxNSV+1;ISV++){
                       	       MBResidual[ISV]=0;
                       	       int EqnBlkCounter = 0;
							   for(J=1;J<ISV-1+1;J++){// -1+1 silliness to agree with VB code
                       	       	       EqnBlkCounter += NEQNsm[J];
                       	       }
                       	       for(J=1;J<NT+1;J++){
                       	       	      int Eqn = EqnBlkCounter+SectnToEqn[ISV][J];
									  Count +=1;
                       	       	      if(SVComp[ISV][J] == 1){
                       	       	      MBResidual[ISV] += c[Eqn]*Compartment[J][7]*Compartment[J][9]*R[Eqn];
                       	       	      }
                       	       }
                       }
                       
 			try {
				File file = new File("MassBalance.csv");
				if(file.exists()){
					file.delete(); // remove residual file
                     		}				
				FileWriter fw = new FileWriter(file.getName(),true);
				BufferedWriter bw = new BufferedWriter(fw); 
				bw.write("State Variable, Initial, In, Transport Out, Volume-Based Sinks Out, Area-Based Sinks Out, Residual, Balance Error ");
				bw.write("\r\n");//line feed
				double ErrorSV = 0;
				double ErrorAcrossSV = 0;
				double SumMBInitial = 0;
				double SumMassIn = 0;
				double SumMassOut = 0;
				for(int ISV = 1; ISV < MaxNSV+1;ISV++){
					SumMBInitial += MBInitial[ISV];
					SumMassIn += SumMassInSV[ISV];
					for(J=1;J<3+1;J++){
						SumMassOut += SumMassOutSV[ISV][J];
					}
					ErrorSV = MBInitial[ISV]+SumMassInSV[ISV]-SumMassOutSV[ISV][1]-SumMassOutSV[ISV][2]-SumMassOutSV[ISV][3]-MBResidual[ISV];
					bw.write(ISV+","+MBInitial[ISV]+","+SumMassInSV[ISV]+","+SumMassOutSV[ISV][1]+","+SumMassOutSV[ISV][2]+","+SumMassOutSV[ISV][3]+","+MBResidual[ISV]+","+ErrorSV);
					ErrorAcrossSV += ErrorSV;
					bw.write("\r\n");//line feed
				}
				
				if(SumMBInitial+SumMassIn == 0){
					System.out.println("Zero initial mass and zero mass in during simulation ...?");
					bw.write("Zero initial mass and zero mass in during simulation ...?");
				}else{
					bw.write("Mass balance error across state variables as a percentage of initial and total input mass = "+ (ErrorAcrossSV * 100) / (SumMBInitial + SumMassIn)+ " %");
				}
				   bw.write("\r\n");//line feed
				bw.close();
			} catch (IOException e){
				 e.printStackTrace();
			}	

                }//end if(DoMassBalance == true)
              
        } // end TimeLoop method
        
        static void ParmUpdate(double Time) {
		
	if (VerboseLevel > 1){
		System.out.println("I'm in ParmUpdate() and Time = "+Time);
	}

	 int I, J, J1, J2, JJ, K, EqnCount, SV1, SV2, Comp1, Comp2, counter;
	 String xstring;
	 float Rate;

        if (!EquationSolver) {
            if (Time==0 || ShellCompartment == true) { // following stuff does not change dynamically
            	    if(ShellCompartment == true){
			WriteTime("Compartment", Time);
			ShellFunction("Compartment.exe");
			}
                  ReadCompartments();
                 
                // EQNTOSECTN(ISV,I), GIVES THE Compartment NUMBER CORRESPONDING TO MASS BALANCE EQUATION I (IN THE Asm MATRIX)for SV ISV
                // SectnToEqn(ISV,I) gives the equation number for compartment i for SV ISV

                for (int i=1; i<MaxNSV+1; i++) {
                	if(VerboseLevel == 3){
                		System.out.println("i = " + i);
                	}
                    counter=0;
                    for (int j=1; j<NT+1; j++) {
                    	    if(VerboseLevel== 3){
                    	    	    System.out.println("j = " + j);
                    	    }
                        SectnToEqn[i][j]= (int) -999; // default
                            EqnToSectn[i][j]= (int) -999; // default
                            if(VerboseLevel == 3){
                            	    System.out.println("SVComp["+i+"]["+j+"] = "+SVComp[i][j]);
                            }
                            if (SVComp[i][j]==1) { // SV I is relevant
                                if (Compartment[j][3]!=0) {
                                	System.out.println("  Error: Compartment[j][3]!=0 in ParmUpdate");
                                	System.exit(1);
                                }
                                    counter++;
                                    EqnToSectn[i][counter] = j; // maps eqn I onto Compartment J
                                    SectnToEqn[i][j] = counter; // maps compartment J onto eqn I
                                    if(VerboseLevel == 3){
                                    	    System.out.println("EqnToSectn["+i+"]["+counter+"] = "+ j);   
                                    	    System.out.println("SectnToEqn["+i+"]["+j+"] = "+ counter);
                                    }
                            }
                    }// end j loop
                }// end i loop
            } // end if Time == 0 || ShellCompartment == true 
            
	       if (Time==0 || ShellLinearKdandTempCoef == true) { 
		       if(ShellLinearKdandTempCoef == true){
			       WriteTime("LinearKdandTempCoef", Time);
			       ShellFunction("LinearKdandTempCoef.exe");
		       }
	       }
                ReadKdandTempCoef();

                // Load V and R arrays
                EqnCount=0;
                counter=0;
                for (int i=1; i<MaxNSV+1; i++) {
                    for (int j=1; j<NT+1; j++) {
                        counter++;
                        if (Compartment[j][3]==0) {
                            // interior compartment
                            if (SVComp[i][j]==1) {
                                // sv I relevant
                                EqnCount++;
                                V[EqnCount]=Compartment[j][7];
                                R[EqnCount] = 1 + Compartment[j][10]*KdandTempCoef[counter][3]/Compartment[j][9]; // 1 + BD*Kd/water content
                            }
                        }
                    }
                }
        
		    if (EqnCount!=NEQN) {
			System.out.println("  Error: EqnCount!=NEQN in ParmUpdate");
						System.exit(1);
					}
		    if (counter != MaxNSVxNT) {
			System.out.println("  Error: counter != MaxNSVxNT in ParmUpdate");
						System.exit(1);
					}                 
	  
		if (Time==0 || ShellInterfaces == true) { 
		       if(ShellInterfaces == true){
			       WriteTime("Interfaces", Time);
			       ShellFunction("Interfaces.exe");
		       }
	       }	
	       ReadInterfaces();

		if (Time==0 || ShellECoefficients == true) { 
		       if(ShellECoefficients == true){
			       WriteTime("Ecoefficients", Time);
			       ShellFunction("Ecoefficients.exe");
		       }
	       }	
	       ReadECoef();       

    		if (Time==0 || ShellFlows == true) { 
		       if(ShellFlows == true){
			       WriteTime("Flows", Time);
			       ShellFunction("Flows.exe");
		       }
	       }	
	       ReadFlows();
    
		if (Time==0 || ShellFlowandEMultipliers == true) { 
		       if(ShellFlowandEMultipliers == true){
			       WriteTime("FlowandEMultipliers", Time);
			       ShellFunction("FlowandEMultipliers.exe");
		       }
	       }	
	       //Note: FlowandEMultipliers.csv is read in Transport() method
		
		if (Time==0 || ShellBoundary == true) { 
		       if(ShellBoundary == true){
			       WriteTime("Boundary", Time);
			       ShellFunction("Boundary.exe");
		       }
	       }	
	       ReadBoundary();
		    
		    if(VerboseLevel == 3){
		    	    System.out.println ("alpha = " +Alpha);
		    }
		    if (CallWiggle == true) {
			UpdateAlpha();
		    }  
            

    } else {
        if (Time==0 || ShellAEqnSolver==true) {
            if (ShellAEqnSolver== true) {
                WriteTime("AEqnSolver",Time);
                ShellFunction("AEqnSolver.exe");
            }
             try {
		       String strFile = externalAfile;
		       BufferedReader br = new BufferedReader( new FileReader(strFile));//create BufferedReader to read csv file
		       String strLine = "";
		       StringTokenizer st = null;
			 for (int i=1; i <NEQN+1;i++) {
				 strLine = br.readLine();
				 st = new StringTokenizer(strLine, ",");
				 for (int j = 1; j < NEQN+1; j++){
					  A[i][j] = Float.parseFloat(st.nextToken());
					 if(VerboseLevel == 3) System.out.println("A["+i+"]["+j+"] = " +A[i][j]);
				 } // end for j
			 } // end for i
			 br.close(); // close input stream
		 }// end try
		 catch(Exception e)
		 {
			System.out.println("Exception while reading "+e);
			System.exit(1);
		 }// end catch
        }//end if (Time==0 || ShellAEqnSolver==true)
         if (Time==0 || GEM.ShellRVEqnSolver==true) {
            if (ShellRVEqnSolver== true) {
                WriteTime("RVEqnSolver",Time);
                ShellFunction("RVEqnSolver.exe");
            }
            if(Dynamic == true){// don't need if in steady state
             	    try {
		       String strFile = externalRVfile;
		       BufferedReader br = new BufferedReader( new FileReader(strFile));//create BufferedReader to read csv file
		       String strLine = "";
		       StringTokenizer st = null;
			 for (int i=1; i <NEQN+1;i++) {
				 strLine = br.readLine();
				 RVEqnSolver[i] = Double.parseDouble(strLine);
				 if(VerboseLevel == 3) System.out.println("RVEqnSolver["+i+"] = " +RVEqnSolver[i]);
			 } // end for i
			 br.close(); // close input stream
		 }// end try
		 catch(Exception e)
		 {
			System.out.println("Exception while reading "+e);
			System.exit(1);
		 }// end catch
            }
         }//end if (Time==0 || GEM.ShellRVEqnSolver==true)
    }// end if (!EquationSolver) 
}// end ParmUpdate() method

static void ReadCompartments() {
	if(VerboseLevel > 1){
		System.out.println("I'm in ReadCompartments() ");
	}

	 try {
		 
		 String strFile = "Compartments.csv";
		 BufferedReader br = new BufferedReader( new FileReader(strFile));//create BufferedReader to read csv file
		 String strLine = "";
		 StringTokenizer st = null;
		 strLine = br.readLine();//"Section,No Adj Sections,Boundary,X (m),Y (m),Z (m),Volume (m^3),Temp (deg C),Water Content (L3w/L3total),Bulk Density (gm/cu m-total)
		 if(VerboseLevel == 3){
		 	 System.out.println("strLine = " + strLine);
		 }
		 for (int i=1; i < NT+1;i++) {
		 	 if(VerboseLevel == 3){
		 	 	 System.out.println("i = " + i);
		 	 }
			 strLine = br.readLine();
			 st = new StringTokenizer(strLine, ",");
			 for (int j = 1; j < NColCompartment+1; j++){
			 	 if(VerboseLevel == 3) System.out.println("j = " + j);
				 Compartment[i][j] = Float.parseFloat(st.nextToken());
				 if(VerboseLevel == 3) System.out.println("Compartment = " +Compartment[i][j]);
			 } // end for j
		 } // end for i
		 br.close(); // close input stream

	 }// end try
	 catch(Exception e)
	 {
		System.out.println("Exception while reading Compartments.csv file: " + e);
	 }// end catch

			 
			 
} // end ReadCompartments() method

static void ReadKdandTempCoef(){
	if(VerboseLevel > 1) System.out.println("I'm in ReadKdandTempCoef() ");
	
	try {
                         	 String strFile = "linearKdandTempCoef.csv";
                         	 BufferedReader br = new BufferedReader( new FileReader(strFile));//create BufferedReader to read csv file
                         	 String strLine = "";
                         	 StringTokenizer st = null;
                         	 strLine = br.readLine();//State Var,Compartment,Linear Kd (cu m-water/gm),Temp Coeff
                         	 if(VerboseLevel == 3) System.out.println("strLine = " + strLine);
                         	 for (int i=1; i < MaxNSVxNT+1;i++) {
                         	 	 if(VerboseLevel == 3) System.out.println("i = " + i);
                         	 	 strLine = br.readLine();
                          	 	 st = new StringTokenizer(strLine, ",");
                         	 	 for (int j = 1; j < NColKdandTempCoef+1; j++){
                         	 	 	 if(VerboseLevel == 3)System.out.println("j = " + j);
                         	 	 	 KdandTempCoef[i][j] = Float.parseFloat(st.nextToken());
                         	 	 	 if(VerboseLevel == 3)System.out.println("KdandTempCoef = " +KdandTempCoef[i][j]);
                         	 	 } // end for j
                         	 } // end for i
                         	 br.close(); // close input stream

			 }// end try
			 catch(Exception e)
			 {
			 	if(VerboseLevel >1) System.out.println("Exception while reading linearKdandTempCoef.csv file: " + e);
			 }// end catch
} // end ReadKdandTempCoef() method

static void ReadInterfaces(){
	if(VerboseLevel>1) System.out.println("I'm in ReadInterfaces() ");
	/*
	' Inface(I,1) = Compartment #
	' REST OF COLUMNS ARE IN GROUPS OF NColPerInface
	' WITH 1 GROUP PER ADJACENT Compartment
	' WITHIN EACH GROUP:
	'     1ST COLUMN IS ADJACENT Compartment #
	'     2ND COLUMN IS INTERFACE AREA IN sq m
	'     3RD COLUMN IS DISTANCE FROM CENTROID OF Compartment TO INTERFACE WITH ADJ Compartment (M)
	'     4TH COLUMN IS MIXING DISTANCE LENGTH (m) TO USE IN LIEU OF DISTANCE BETWEEN CENTROIDS
	'     5th column is alpha
	*/
	int Count;
	try {
	       String strFile = "Interfaces.csv";
	       BufferedReader br = new BufferedReader( new FileReader(strFile));//create BufferedReader to read csv file
	       String strLine = "";
	       StringTokenizer st = null;
	       strLine = br.readLine();//Section,Adj Section,Interface Area (m^2),Section Distance to Interface of Adj Section (m),Mixing Length (m),Adj Section,Interface Area (m^2),Section Distance to Interface of Adj Section (m),Mixing Length (m),Adj Section,Interface Area (m^2)
	       if(VerboseLevel == 3) System.out.println("strLine = " + strLine);
	       for (int i=1; i < NT+1;i++) {
                         	 	 if(VerboseLevel == 3) System.out.println("i = " + i);
                         	 	 // J1 WILL BE # OF CompartmentS ADJ TO Compartment I
                         	 	 int j1 = (int) Compartment[i][2];
                         	 	 // J2 WILL BE # OF COLUMNS TO READ IN FIELD I (Compartment I)
                         	 	 int j2 = j1 * NColPerInface + 1;
                         	 	 strLine = br.readLine();
                          	 	 st = new StringTokenizer(strLine, ",");
                          	 	 Count = 0;
                         	 	 for (int j = 1; j < j2+1; j++){
                         	 	 	 if(j == 2) Count = 0;
                         	 	 	 Count += 1;
                         	 	 	 if(VerboseLevel == 3) System.out.println("j = " + j);
                         	 	 	 Inface[i][j] = Float.parseFloat(st.nextToken());
                         	 	 	 if(j==1 && Inface[i][j] != i) {
							System.out.println("data error in Interfaces.csv"); 
							System.exit(1);
                         	 	 	 }
                         	 	 	 if(VerboseLevel == 3) System.out.println("Inface = " +Inface[i][j]);
                         	 	 	 if(Count == NColPerInface){
                         	 	 	 	 if (Inface[i][j] != 1) CallWiggle = true;
                         	 	 	 }
                         	 	 	 if(Count == NColPerInface)Count = 0;
                         	 	 } // end for j
                         	 } // end for i
                         	 br.close(); // close input stream
			 }// end try
			 catch(Exception e)
			 {
			 	System.out.println("Exception while reading Interfaces.csv file: " + e);
			 }// end catch
} // end ReadInterfaces() method

static void ReadECoef(){
	/*
	'Read dispersion coefficients
	' ECoef(I,1) = State Variable
	' ECoef(I,2) = Compartment I
	' REST OF COLUMNS ARE IN GROUPS OF NColPerEcoef
	' WITH 1 GROUP PER ADJACENT Compartment
	' WITHIN EACH GROUP:
	'     1ST COLUMN IS ADJACENT Compartment #
	'     2nd COLUMN IS Dispersion Coefficient in sq m/day
	*/
	if(VerboseLevel > 1) System.out.println("I'm in ReadECoef() ");
	try {
		       String strFile = "ECoefficients.csv";
		       BufferedReader br = new BufferedReader( new FileReader(strFile));//create BufferedReader to read csv file
		       String strLine = "";
		       StringTokenizer st = null;
		       strLine = br.readLine();//State Variable,Compartment,Adj Compartment,E (m^2/day),Adj Compartment,E (m^2/day)
		       if(VerboseLevel == 3) System.out.println("strLine = " + strLine);
		       int counter = 0;
		       for (int i=1; i < MaxNSV+1;i++) {
				 if(VerboseLevel == 3) System.out.println("i = " + i);
				 for (int j = 1; j < NT+1; j++){
					  strLine = br.readLine();
					  st = new StringTokenizer(strLine, ",");
					 if(VerboseLevel == 3) System.out.println("j = " + j);
					 counter++;
					 int j1 = (int) Compartment[j][2];// J1 WILL BE # OF CompartmentS ADJ TO Compartment I
					 int j2 = j1 * GEM.NColPerEcoef + 2;// J2 WILL BE # OF COLUMNS TO READ IN record counter
					 for (int k = 1; k < j2+1; k++){
						 if(VerboseLevel == 3) System.out.println("k = " + k);
						ECoef[counter][k] = Float.parseFloat(st.nextToken());
						if(VerboseLevel == 3) System.out.println("ECoef = " +ECoef[counter][k]); 	 
					 }
				 } // end for j
			 } // end for i
			 br.close(); // close input stream
		 }// end try
				 catch(Exception e)
		 {
			System.out.println("Exception while reading ECoefficients.csv file: " + e);
		 }// end catch
	
		} // end ReadECoef() method
	
	static void ReadFlows(){
	/*
		'Flow(I,1) = state variable
	    ' Flow(I,2) = Compartment I
	     ' REST OF COLUMNS ARE IN GROUPS OF NColPerFlow
	     ' WITH 1 GROUP PER ADJACENT Compartment
	     ' WITHIN EACH GROUP:
	     '     1ST COLUMN IS ADJACENT Compartment #
	     '     2nd COLUMN IS flow in cu m/day
	 */
	if(VerboseLevel > 1) System.out.println("I'm in ReadFlows() ");
	try {
		       String strFile = "Flows.csv";
		       BufferedReader br = new BufferedReader( new FileReader(strFile));//create BufferedReader to read csv file
		       String strLine = "";
		       StringTokenizer st = null;
		       strLine = br.readLine();//Compartment,Adj Compartment,Flow (cmd),Adj Compartment,Flow (cmd)
		       if(VerboseLevel == 3) System.out.println("strLine = " + strLine);
		       int RowCount = 0;
		       for(int ISV = 1;ISV < MaxNSV+1;ISV++){
			       for (int i=1; i < NT+1;i++) {
			       	       RowCount += 1;
					 if(VerboseLevel == 3) System.out.println("i = " + i);
					 // J1 WILL BE # OF CompartmentS ADJ TO Compartment I
					 int j1 = (int) Compartment[i][2];
					 // J2 WILL BE # OF COLUMNS TO READ IN FIELD I (Compartment I)
					 int j2 = j1 * GEM.NColPerFlow + 2;
					 strLine = br.readLine();
					 st = new StringTokenizer(strLine, ",");
					 for (int j = 1; j < j2+1; j++){
						 if(VerboseLevel == 3) System.out.println("j = " + j);
						 Flow[RowCount][j] = Float.parseFloat(st.nextToken());
						 if (j == 2 && Flow[RowCount][j]!= i){
							System.out.println("  Error: j == 1 && Flow[i][j]!= i) in ReadFlows()");
							System.exit(1);
						 }
						if(VerboseLevel == 3 )System.out.println("Flow = " +Flow[RowCount][j]);
					 } // end for j
				 } // end for i
		       }//end for ISV
			 br.close(); // close input stream
		 }// end try
				 catch(Exception e)
		 {
			System.out.println("Exception while reading Flows.csv file: " + e);
		 }// end catch
	
		} // end ReadFlow() method	

	 static void ReadBoundary(){
 	 
	if(VerboseLevel> 1) System.out.println("I'm in ReadBoundary() ");
	try {
		       String strFile = "Boundary.csv";
		       BufferedReader br = new BufferedReader( new FileReader(strFile));//create BufferedReader to read csv file
		       String strLine = "";
		       StringTokenizer st = null;
		       strLine = br.readLine();//Compartment,Adj Compartment,Flow (cmd),Adj Compartment,Flow (cmd)
		       if(VerboseLevel == 3) System.out.println("strLine = " + strLine);
		       int counter = 0;
		       for (int i=1; i < MaxNSVxNT+1;i++) {
				 if(VerboseLevel == 3) System.out.println("i = " + i);
				 strLine = br.readLine();
				 st = new StringTokenizer(strLine, ",");
				 for (int j = 1; j < NColBoundary+1; j++){
					 if(VerboseLevel == 3) System.out.println("j = " + j);
					 Boundary[i][j] = Float.parseFloat(st.nextToken());
					 if(VerboseLevel == 3) System.out.println("Boundary = " +Boundary[i][j]);
				 } // end for j
			 } // end for i
			 br.close(); // close input stream
		 }// end try
				 catch(Exception e)
		 {
			System.out.println("Exception while reading Boundary.csv file: " + e);
		 }// end catch

	} // end ReadBoundary() method
	
	static void UpdateAlpha(){
	if(VerboseLevel> 1) System.out.println("I'm in UpdateAlpha() ");
	
	int I,IC1,IC2,J,K,L,NumAdjSectn1,NumAdjSectn2,AdjSectnNum;
	float Dist1,Dist2;
	
	for(I=1;I<NT+1;I++){
		AlphaWeighted[I][1] = (float) I;
		if (Compartment[I][3] == 0) {// I is an interior compartment
			if(VerboseLevel == 3) System.out.println("Interior compartment "+I);
			NumAdjSectn1 = (int) Compartment[I][2];
			if(VerboseLevel == 3) System.out.println("NumAdjSectn1 = "+NumAdjSectn1);
			IC1 = 2;//initialize
			IC2 = 2;//initialize
			for(K = 1;K<NumAdjSectn1+1;K++){//walk adjacent compartments to compartment I
				if(K == 1){
					IC1 = 2;
				}else{
					IC1 = IC1+NColPerInface;
				}
				AdjSectnNum = (int)Inface[I][IC1];
				if(VerboseLevel == 3) System.out.println(" Adj Comp ="+AdjSectnNum);
				if(AdjSectnNum == I) {
					System.out.println("  Error in UpdateAlpha, AdjSectnNum = I");
					System.exit(1);	
				}
				Dist1 = Inface[I][IC1 + 2]; //distance from centroid of I to interface with AdjSectnNum
				if(Dist1 == -999) {
					System.out.println("  Error in UpdateAlpha, Dist1 = -999");
					System.exit(1);	
				}
				//now get distance from AdjSectnNum centroid to interface of I,AdjSectnNum
				NumAdjSectn2 = (int)Compartment[AdjSectnNum][2];
				if(VerboseLevel == 3) System.out.println("NumAdjSectn2 "+NumAdjSectn2);
				Dist2 = -999;//initialize
				for(L = 1;L<NumAdjSectn2+1;L++){
					if(L == 1){
						IC2 = 2;
					}else{
						IC2 = IC2+NColPerInface;
					}
					if((int)Inface[AdjSectnNum][IC2] == I) { 
						Dist2 = Inface[AdjSectnNum][IC2 + 2];
						if(Dist2 == -999) {
							System.out.println("  Error in UpdateAlpha, Dist2 = -999");
							System.exit(1);	
						}
                     			}
				}//end L loop
				AlphaWeighted[I][2 * K] = AdjSectnNum;
				AlphaWeighted[I][2 * K + 1] = Dist2 / (Dist1 + Dist2);
				if(VerboseLevel == 3) System.out.println("AlphaWeighted = "+AlphaWeighted[I][2 * K + 1]);
				if(AlphaWeighted[I][2*K+1] > 1 || AlphaWeighted[I][2*K+1] < 0){
					System.out.println("  Error in UpdateAlpha, AlphaWeighted > 1 or < 0)");
					System.exit(1);
				}
			}//end K loop
		}//end if (Compartment[I][3] == 0
	}// end I loop
} // end UpdateAlpha() method

	 
	static void LoadATrans() {
	
	// Method TO SET UP MASS BALANCE EQUATIONS IN COEFFICIENT MATRIX A AND LHS VECTOR BBC
	int I,J,IStart, JStart, IA, JA, ISV;
	if(VerboseLevel > 1) System.out.println("I'm in LoadATrans NEQN = " + NEQN);
	
	// initialize arrays
	for (I = 1; I < NEQN+1; I++){
		BBC[I] =  0.0;
		for (J = 1; J < NEQN+1; J++){
			A[I][J] =  0.0;
		} // end for J
	}// end for I
	
	//MB code!
	if(DoMassBalance == true){
		MBBoundaryRowCount = 0;//initialize
    	}

	// PUT TRANSPORT TERMS INTO A AND BBC
	
	for (ISV = 1; ISV < MaxNSV+1; ISV++) {	
		// first calculate number of equations for SV ISV
		NEQNsm[ISV] = 0;
		for (J = 1; J < NT+1; J++) {
			if (Compartment[J][3] == 0) {// interior compartment
				NEQNsm[ISV] +=  SVComp[ISV][J];
			} // end if
		} // end for J
		
		if (NEQNsm[ISV] == 0) {
			System.out.println("  Error In LoadATrans(): NEQNsm[ISV] == 0");
			System.exit(1);
		}
		if(VerboseLevel == 3) System.out.println("I'm in LoadATrans NEQNsm[" + ISV + "] = " + NEQNsm[ISV]);
		// Asm and Bsm WILL BE SUBSETS OF THE LARGER SIMULTANEOUS SYSTEM
		Asm = new float[NEQNsm[ISV]+1][NEQNsm[ISV]+1];	
		Bsm = new float[NEQNsm[ISV]+1];
		// initialize
		for (I = 1; I < NEQNsm[ISV]+1; I++){
			Bsm[I] = (float) 0.0;
			for (J = 1; J < NEQNsm[ISV]+1; J++){
				Asm[I][J] = (float) 0.0;
			} // end for J
		}// end for I
			
		Transport(ISV);// call Transport method in this Class
		
		IStart = 1; //initialize for ISV = 1
		JStart = 1; // initialize for ISV = 1
		for (I=1;I<=ISV;I++){
				if(I>1){
				IStart += NEQNsm[I-1];
				JStart = IStart;
			}
		}
	
		//if (ISV == 1) {
			//IStart = 1;
			//JStart = 1;
		//} else {
			//IStart += NEQNsm[ISV-1];
			//JStart = IStart;
		//} // end else
		
	
		IA = 0; // Row counter of matrix Asm
		for (I = IStart; I < IStart + NEQNsm[ISV]-1+1; I++) {// Keith note: "-1+1" foolishness to maintain comparison with VB code
			IA += 1;
			JA = 0; // column counter for matrix Asm
			//System.out.println("I'm in LoadATrans and I = " + I);//here!
			for (J = JStart; J < JStart + NEQNsm[ISV]-1+1; J++) {
				JA += 1;
				A[I][J] = Asm[IA][JA];
			} //end for J
			BBC[I] = Bsm[IA];
		} // end for I
	} //end for ISV


} // end LoadATrans() method

 static void LoadALinearSrcSnk() {
	
	int ARow, ACol, Comp, Comp1, Comp2, SV, SV1, SV2, EqnBlkCounter,I,J,KdandTempCoefCounter;
	float Area, Temp, Rate;
	double Dummy;

	//12_14_16 change:
	float WCcomp1, WCint; //water content compartment 1 and water content of comp1,comp2 interface
	
	if(VerboseLevel > 1) System.out.println("I'm in LoadALinearSrcSnk method");
	
	// MB code!
	if(DoMassBalance == true){
	    MBExtVolSinksRowCount = 0;//initialize
	    MBExtAreaSinksRowCount = 0;//initialize
	    for(I =1;I<MaxNSVxNT+1;I++){
	    	    for(J=1;J<3+1;J++){
	    	    	 MBExtVolSinks[I][J] = -999;//initialize
	    	    	 MBExtAreaSinks[I][J] = -999;//initialize
	    	    	 if(J == 3) MBExtAreaSinks[I][J+1] = -999;//initialize -- get 4th column
	    	    }
	    }
	}
	// read VolumeSrcSnks.csv
	try{
		String strFile = "VolumeSrcSnks.csv";
		BufferedReader br = new BufferedReader( new FileReader(strFile));//create BufferedReader to read csv file
		String strLine = "";
		StringTokenizer st = null;
		strLine = br.readLine();//Compartment,First SV,Second SV,Rate Constant (/day) wrt First SV
		if(VerboseLevel == 3) System.out.println(strLine);
	
		do {
			strLine = br.readLine();//data record
			if(VerboseLevel == 3) System.out.println(strLine);
			st = new StringTokenizer(strLine, ",");
			Comp = Integer.parseInt(st.nextToken());
			if (Comp == -999) break;// end of file
			SV1 = Integer.parseInt(st.nextToken());
			SV2 = Integer.parseInt(st.nextToken());
			Rate = Float.parseFloat(st.nextToken());
			if (Compartment[Comp][3] == 0){
				
				EqnBlkCounter = 0; //initialize -- counts rows of A matrix
				KdandTempCoefCounter = 0;//initialize -- counts rows of KdandTempCoef array
				for (J = 1; J < SV1-1+1; J++) {// Keith note: "-1+1" foolishness to maintain comparison with VB code
					EqnBlkCounter += NEQNsm[J];
					KdandTempCoefCounter += NT;
				}
				ARow=EqnBlkCounter+SectnToEqn[SV1][Comp];
				KdandTempCoefCounter += Comp;

				if (SV2 == 1) {
					if(VerboseLevel == 3) System.out.println("SV2= "+SV2+"Comp = "+Comp);
					if(Rate < 0) {//sink goes on main diagonal
						ACol = ARow;
					} else{
						ACol = SectnToEqn[SV2][Comp];
						}
				}else{
					if(Rate<0) {//sink goes on main diagonal
						ACol = ARow;
					}else{
						ACol = NEQNsm[SV2-1]+SectnToEqn[SV2][Comp];
					}
				} // end if (SV2 == 1)	
				Temp = Compartment[Comp][8];
				//following line assumes that only dissolved conc is undergoing rate kinetics
				Dummy = V[ARow] * Rate * Math.pow(KdandTempCoef[KdandTempCoefCounter][4],(Temp - 20));
				//following line assumes that the total conc (dissolved + sorbed) is undergoing rate kinetics
				// Dummy = R(ARow) * V(ARow) * Rate * KdandTempCoef(KdandTempCoefCounter, 4) ^ (Temp - 20);
				
				A[ARow][ACol] += Dummy;
				//System.out.println("after s/s change A["+ARow+"]["+ACol+"] = "+A[ARow][ACol]);
				
				//MB code
				if(DoMassBalance == true){
					if(SV1 == SV2){//external source/sink, i.e. not an across-SV transfer
						MBExtVolSinksRowCount += 1;
						if(MBExtVolSinksRowCount > MaxNSVxNT){
							System.out.println("Error: number of records in mass balance array MBExtVolSinks exceeds dimensioned length of array");
							System.exit(1);
						}
						MBExtVolSinks[MBExtVolSinksRowCount][1] = SV1;
						MBExtVolSinks[MBExtVolSinksRowCount][2] = Comp;
						MBExtVolSinks[MBExtVolSinksRowCount][3] = Dummy * Compartment[Comp][9];						
					}
				}//end if(DoMassBalance == true)
			} else{
				System.out.println("In VolumeSrcSnks.csv a sink is specified for a boundary compartment -- being ignored");
			}	
		} while (Comp != -999  );
		br.close(); // close input stream
	} // end try
	catch(Exception e)
	{
		System.out.println("Exception while reading VolumeSrcSnks.csv file: " + e); 
	} // end catch
	
	// read AreaSrcSnks.csv
	try{
		String strFile = "AreaSrcSnks.csv";
		BufferedReader br = new BufferedReader( new FileReader(strFile));//create BufferedReader to read csv file
		String strLine = "";
		StringTokenizer st = null;
		strLine = br.readLine();//SV,First Compartment, Second Compartment, Rate Const (/day) wrt First Compartment
		if(VerboseLevel == 3) System.out.println(strLine);
		do {
			strLine = br.readLine();//data record
			st = new StringTokenizer(strLine, ",");
			SV = Integer.parseInt(st.nextToken());
			if(VerboseLevel == 3) System.out.println("SV = "+SV);
			if (SV == -999) break;
			Comp1 = Integer.parseInt(st.nextToken());
			Comp2 = Integer.parseInt(st.nextToken());
			Rate = Float.parseFloat(st.nextToken());
			if(VerboseLevel == 3) System.out.println("Comp1, Comp2, Rate = "+Comp1+","+Comp2+","+Rate);
			EqnBlkCounter = 0;//initialize
			for (J = 1; J < SV-1+1; J++) {// Keith note: "-1+1" foolishness to maintain comparison with VB code
				EqnBlkCounter += NEQNsm[J];
			}
			ARow = EqnBlkCounter+SectnToEqn[SV][Comp1];
			
			ACol = ARow; //initialize
			if (Rate<=0){//sink goes on main diagona'
				ACol = ARow;
			}else{
				if (Compartment[Comp2][3]!=0){
					System.out.println("Warning: GEM not set up to accept an area-based source from a boundary/dummy compartment to an interior compartment user should specify as external loading instead in Loads.csv"); 	
				}else{
					ACol = EqnBlkCounter+SectnToEqn[SV][Comp2];
				}
			}
			Area = -999; // initialize
			Area = GetArea(Comp2, Comp1); // get interface area
			if(VerboseLevel == 3) System.out.println("Area = "+Area);
			if(VerboseLevel == 3) System.out.println("Arow = "+ARow+" ACol = "+ACol);
			if (Area == -999) {
				System.out.println("Error!  interface Area = -999");
				System.exit(1);
			}

			//12_14_16 change:
			Dummy = -999;//initialize
            	if (Compartment[Comp2][3] != 0){
 				//Comp2 is boundary or dummy compartment
                		Dummy = Area * Rate;
            	}else{
				//Comp2 is interior compartment and we need interface water content
                		WCcomp1 = Compartment[Comp1][9];
                		WCint = WCcomp1;
                		if (Compartment[Comp2][9] < WCint) { 
                    		WCint = Compartment[Comp2][9]; // minimum water content of comp1,comp2
					if (WCint == -999){
						System.out.println("Error!  Water content at Comp1, Comp2 interface is -999 in LoadALinearSrcSnk");
						System.exit(1);
					}
				}
				Dummy = Area * Rate * WCint / WCcomp1;
                	}
			//Dummy = Area*Rate;
			
			if (Dummy == -999) {
				System.out.println("Error!  Dummy = -999 in LoadALinearSrcSnk");
				System.exit(1);
			}
			
			A[ARow][ACol] += Dummy;
			//MB code!
			if(DoMassBalance == true){
				if(Compartment[Comp2][3] !=0 && Rate<0){//external sink to comp2 a boundary/dummy compartment
					MBExtAreaSinksRowCount+= 1;
					if(MBExtAreaSinksRowCount > MaxNSVxNT){
							System.out.println("Error: number of records in mass balance array MBExtAreaSinks exceeds dimensioned length of array");
							System.exit(1);
					}
					MBExtAreaSinks[MBExtAreaSinksRowCount][1] = SV;
					MBExtAreaSinks[MBExtAreaSinksRowCount][2] = Comp1;
					MBExtAreaSinks[MBExtAreaSinksRowCount][3] = Comp2;
					MBExtAreaSinks[MBExtAreaSinksRowCount][4] = Dummy * Compartment[Comp1][9];
                		}
            		}

		} while (SV != -999  );
		br.close(); // close input stream
	} // end try
	catch(Exception e)
	{
		System.out.println("Exception while reading AreaSrcSnks.csv file: " + e); 
	} // end catch

}// end LoadALinearSrcSnk() method

	static void Aout(float TempTime){//KWL 3/15 change
		if(VerboseLevel >1) System.out.println("I'm in Aout() ");
		try {
			FileWriter fw = new FileWriter("Aout.dng");
			BufferedWriter out = new BufferedWriter(fw);
			out.write("At time "+TempTime);
			out.write("\r\n");// line feed
			for (int I = 1; I < NEQN+1;I++){
				for (int J = 1; J < NEQN+1;J++){
					String dummy = Double.toString(A[I][J]);
 					out.write(dummy);
 					if (J < NEQN){
 						out.write(",");
 					}
				}
			out.write("\r\n");// line feed
			}
			out.close();
		} catch (IOException e){
			 e.printStackTrace();
		}
	} // end Aout method
	
	static void Bout(float TempTime){//KWL 3/15 change
		if(VerboseLevel> 1) System.out.println("I'm in Bout() ");
		try {
			FileWriter fw = new FileWriter("Bout.dng");
			BufferedWriter out = new BufferedWriter(fw);
			out.write("At time "+TempTime);
			out.write("\r\n");// line feed
			for (int I = 1; I < NEQN+1;I++){
				String dummy = Double.toString(b[I]);
 					out.write(dummy);
 			out.write("\r\n");// line feed
			}
			out.close();
		} catch (IOException e){
			 e.printStackTrace();
		}
	} // end Bout method
	
	
	static void RVout(float TempTime){// KWL 3/14 change
		if(VerboseLevel > 1 )System.out.println("I'm in RVout() ");
		try {
			FileWriter fw = new FileWriter("RVout.dng");
			BufferedWriter out = new BufferedWriter(fw);
			out.write("At time "+TempTime);
			out.write("\r\n");// line feed			
			out.write("The diagonal elements of RV matrix are ...");
			out.write("\r\n");// line feed
			for (int I = 1; I < NEQN+1;I++){
				if (EquationSolver == false){
					String dummy = Double.toString(R[I]*V[I]);
 					out.write(dummy);
				}else{
					System.out.println("VB code needs translating here in RVout");
					//Print #1, RVEqnSolver(I)
				}
			out.write("\r\n");// line feed
			}
			out.close();
		} catch (IOException e){
			 e.printStackTrace();
		}
	} // end RVout method
	

	static float GetArea(int I, int J) {
		float Area = (float) -999.0;
		if(VerboseLevel > 1) System.out.println("in GetArea() method ");
		int K, IC;
		if(VerboseLevel == 3) System.out.println("I = "+I+" J = "+J);
		int NAS = (int) Compartment[I][2]; // number of compartments adjacent to compartment I
		if(VerboseLevel == 3) System.out.println("NAS = "+NAS);
		Boolean found = false;
		IC = -999; //initialize
		for (K = 1; K < NAS+1; K++){
			if(VerboseLevel == 3) System.out.println("K = "+K);
			if (K == 1){
				 IC = 2;
			} else {
				 IC += NColPerInface;
				 if(VerboseLevel == 3) System.out.println("NColPerInface = "+NColPerInface);
				 if(VerboseLevel == 3) System.out.println("IC = "+IC);
			}
			if (IC == -999) {
				System.out.println("  Error in GetArea ! IC = -999");
				System.exit(1);
			}
			if(VerboseLevel == 3) System.out.println("Inface["+I+"]["+IC+"] = "+Inface[I][IC]);
			if (Inface[I][IC] == J) {
				Area = Inface[I][IC+1];
				if(VerboseLevel == 3) System.out.println("Area = "+Area);
				found = true;
				if(VerboseLevel == 3) System.out.println("found = true in GetArea");
				break; // exit K loop
			}
		} // end K loop
		if (found == false) {
			System.out.println("  Error in GetArea method -- didn't find interface area for compartments "+I+" and "+J);
			System.exit(1);
		}
		return Area;
	} // end GetArea() method
	
	static void ProfileSnapShot(float TempTime){
		if(VerboseLevel > 1) System.out.println("I'm in ProfileSnapShot() ");
		int IC,ISV,I;
		
		try {
			FileWriter fw = new FileWriter("ProfileSnapShot.csv");
			BufferedWriter out = new BufferedWriter(fw);
			if (EquationSolver == false) {
				out.write("Time,State Variable,Compartment,X,Y,Z,Concentration");
				out.write("\r\n");// line feed	
				IC = 0;
				for (ISV = 1; ISV < MaxNSV+1;ISV++){
					for (I = 1; I < NEQNsm[ISV]+1;I++){
						IC++;
						//start KWL 3/6/12 changes
						int Comp = EqnToSectn[ISV][I];
						out.write(TempTime+","+ISV+","+Comp+","+Compartment[Comp][4]+",");
						if (Linear == false && Explicit == false){
							out.write(Compartment[Comp][5]+","+Compartment[Comp][6]+","+c[IC]/NewtonSF);
						}else{
							out.write(Compartment[Comp][5]+","+Compartment[Comp][6]+","+c[IC]);
						}
						//end KWL 3/6/12 changes
					out.write("\r\n");// line feed	
					}
				}
			} else {// EquationSolver == true
				/*
				Print #1, "TIME"; ","; "VARIABLE"; ","; "VALUE"
				For I = 1 To NEQN
				     If Linear = False And Explicit = False Then
					Print #1, TempTime; ","; I; ","; c(I) / NewtonSF
				     Else
					Print #1, TempTime; ","; I; ","; c(I)
				    End If
				Next I
			    End If	
			    */
			    out.write("Time,Variable,Value");
				out.write("\r\n");// line feed	
				for (I = 1; I < NEQN+1;I++){
					out.write(TempTime+","+I+",");
					if (Linear == false && Explicit == false){
						String stg = Double.toString(c[I]/NewtonSF);
						out.write(stg);
					}else{
						String stg = Double.toString(c[I]);
						out.write(stg);
					}
				out.write("\r\n");// line feed	
				}
			}
			out.close();
		} catch (IOException e){
			 e.printStackTrace();
		}
	}// end ProfileSnapShot method
	
	static void OneCompTimeSeries(int NumSectnPrint, float TempTime){
		if(VerboseLevel > 1) System.out.println("I'm in OneCompTimeSeries() ");
		int IC,I,J;
		    java.io.File outfile;
		
		try {
			File file = new File("OneCompTimeSeries.csv");
			if (TempTime <= DelT){
				if(file.exists()){
					file.delete(); // remove residual file
				}
			}
		   FileWriter fw = new FileWriter(file.getName(),true);
                    BufferedWriter bw = new BufferedWriter(fw); 
                    if (TempTime <= DelT){
                        bw.write("Time, ");
                        for (I = 1;I < MaxNSV+1;I++){
                        	if(SVComp[I][NumSectnPrint]==1){
                        		bw.write("SV("+I+"),");
                            	}
                        }
                           bw.write("\r\n");//line feed
                    }
                    if (EquationSolver == false) {
                    	    bw.write(TempTime+",");
			for (I = 1; I < MaxNSV+1;I++){
				if (SVComp[I][NumSectnPrint]==1){
					IC = 0;
					for (J=1;J<I-1+1;J++){// Keith note: "-1+1" foolishness to maintain comparison with VB code
						IC += NEQNsm[J];
					}
					//start KWL 3/6/12 changes
					int Equation = SectnToEqn[I][NumSectnPrint];
					if (Linear == false && Explicit == false){
						String dummy = Double.toString(c[Equation+IC]/NewtonSF);
						bw.write(dummy+",");
					} else{
						String dummy = Double.toString(c[Equation+IC]);
						bw.write(dummy+",");
					}
					//end KWL 3/6/12 changes
				}// end if (SVComp[I][NumSectnPrint]==1)
			}// I loop
		} else {// EquationSolver == true
			System.out.println("Error: OneCompTimeSeries not coded for EquationSolver = true in ProfileSnapShot");
			System.exit(1);
		}
			bw.write("\r\n");// line feed	
			bw.close();
		} catch (IOException e){
			 e.printStackTrace();
		}
	}// end OneCompTimeSeries method
	
	static void AllTimeSeries(float TempTime){
		if(VerboseLevel > 1) System.out.println("I'm in AllTimeSeries(), TempTime = "+TempTime);
		int I;
                java.io.File outfile;

		try {
                    File file = new File("AllTimeSeries.csv");
                    
                    if (TempTime <= DelT){
                        if(file.exists()){
                        file.delete(); // remove residual file
                        }
                    }
                    FileWriter fw = new FileWriter(file.getName(),true);
                    BufferedWriter bw = new BufferedWriter(fw); 
                    if (TempTime <= DelT){
                        bw.write("Time, ");
                        for (I = 1;I < NEQN+1;I++){
                            bw.write("c("+I+"),");
                        }
                           bw.write("\r\n");//line feed
                    }
                    bw.write(TempTime+",");
                    for (I = 1;I < NEQN+1;I++){
                            if (Linear == false && Explicit == false){
                            	    //start KWL 3/6/12 changes
                            	    String dummy = Double.toString(c[I]/NewtonSF);
                            	    	    bw.write(dummy+",");
                             }else{
                             	     String dummy = Double.toString(c[I]);
                             	     bw.write(dummy+",");        
                            }
                            //end KWL 3/6/12 changes
                    }
                 bw.write("\r\n");
                 bw.close();
		} catch (IOException e){
			 e.printStackTrace();
		}	
		
	}// end AllTimeSeries method
	
	static void Transport(int ISV)  {

	if(VerboseLevel > 1) System.out.println("I'm in Transport() method ");
	int I,J,K;
	int IK, JC, NAS, ICInface, ICEcoef, ICFlow, RowCount, IADJ, IAS, JJ, JStart, NASINT, SV, Comp,Row, Col, CompRow, CompCol;
	float Q, QMult, SumOutflow, SumInflow, T, Tadv, Tdisp,EPrime, AvgAdjL, AvglNew, E, EMult, Avgl, Area, AlphaToUse,WCk;
	float Multiplier;
	
	// 12_14_16 change
	float WCint;// water content at interface of K and IAS
	WCint = -999;//initialize
	
	if(VerboseLevel == 3) System.out.println("Considering State Variable ISV = " +ISV);
	RowCount = (ISV - 1) * NT; //' will count rows of all MaxNSVxNT arrays
	double CiCoef = 0;//initialize
	double CjCoef = 0;//initialize
	for (K = 1; K < NT+1; K++) { // walk compartments
		if(VerboseLevel == 3) System.out.println("Compartment K = " +K);
		RowCount = RowCount + 1;
		if(VerboseLevel == 3) System.out.println("Row count of MaxNSVxNT arrays = " +RowCount);
		if (Compartment[K][3] == 0 && SVComp[ISV][K] == 1) {// Interior compartments only and SV ISV relevant to compartment K
	
			WCk = Compartment[K][9]; // water content compartment K
			//  NOW FIND CompartmentS ADJACENT TO Compartment K
			// NAS WILL BE # OF ADJACENT CompartmentS (INCLUDING BOUNDARY CompartmentS)
			NAS = (int) Compartment[K][2];
			if (NAS ==  -999) {
				System.out.println("  Error: NAS = -999");
				System.exit(1);
			}
	
			// ICInface WILL BE COLUMN COUNTER FOR ARRAY Inface
			// ICTransParms will be column counter for array TransParms
			SumOutflow = 0; // FLOW MASS BALANCE CHECK FOR NON-BOUNDARY CompartmentS
			SumInflow = 0;	
	
			ICInface = 2;// initialize
			ICFlow = 3; // initialize
			ICEcoef = 3; // initialize
			for (J = 1; J < NAS +1; J++) { // walk adjacent compartments
				
				if (J > 1) {
					ICInface += NColPerInface;
					ICFlow += GEM.NColPerFlow;
					ICEcoef += GEM.NColPerEcoef;
				}// end if
				
				// IAS WILL BE ADJACENT Compartment #
				IAS = (int) Inface[K][ICInface];
				if(VerboseLevel == 3) System.out.println("Adj Compartment IAS = " +IAS);
				
				if  (IAS == -999) {
					System.out.println("  Error: IAS = -999");
					System.exit(1);
					}
				if  (IAS == K) {
					System.out.println("  Error: IAS = K");
					System.exit(1);
				}
				if (Flow[RowCount][ICFlow] != IAS){
					System.out.println("Error: Flow[K][ICFlow] != IAS in Transport");
					System.exit(1);
				}
				
				// 12_14_16 change
				if (Compartment[IAS][3] == 0 ){//IAS is interior Compartment
					WCint = WCk;
					if (Compartment[IAS][9] < WCk ){
						WCint = Compartment[IAS][9];
					}
				}

				Q = Flow[RowCount][ICFlow + 1];
				if(VerboseLevel == 3) System.out.println("Q = " +Q);
				if (Q > 0) {
					SumInflow += Q;
				}
				else {
					SumOutflow +=  Math.abs(Q);
				} // end if
				if(Inface[K][ICInface] != IAS){
					System.out.println("Error: Inface[K][ICInface] != IAS in Transport");
					System.exit(1);
				}
				Area = Inface[K][ICInface + 1];
				if(VerboseLevel == 3) System.out.println("Area = " +Area);
				if(ECoef[RowCount][ICEcoef] != IAS){
					System.out.println("Error: ECoef[RowCount][ICEcoef] != IAS in Transport");
					System.exit(1);
				}
				E = ECoef[RowCount][ICEcoef + 1];
				if(VerboseLevel == 3) System.out.println("E = " +E);
				if (E != 0 && E != -999){
					if (MixLengthOption == 1) { // use x,y,z centroid method
						MixLength = GetLength(1,K,IAS); //GET mixing length DISTANCE BETWEEN CompartmentS K AND IAS
						//System.out.println("MixLength = "+MixLength);
					}
					else { // user-provided mixing length override
						MixLength = Inface[K][ICInface + 3];
					} // end if (MixLengthOption == 1)
					if (Area == -999) {
						System.out.println("  Error: Area = -999");
						System.exit(1);
					}
					if (MixLength == -999) {
						System.out.println("  Error: MixLength = -999");
						System.exit(1);
					}
					 EPrime = E*Area/MixLength;
				}
				else { // E not a transport process
					EPrime = 0;
				} // end if (E != 0 && E != -999)
				if(VerboseLevel == 3) System.out.println("EPrime = " +EPrime);
				
				//AlphaToUse = Alpha; // default
				AlphaToUse = Inface[K][ICInface + 4];
				if(VerboseLevel == 3) System.out.println("AlphaToUse = " +AlphaToUse);
	
				// T will be transport terms associated with adjacent compartment IAS
				// First consider off-diagonal elements
				EMult = (float) 1.0; // default
				QMult = (float) 1.0; // default
				if (CallFlowandEMultipliers == true) EMult = ReadFlowandEMultipliers( "EMult",ISV, K, IAS, K, IAS);
				if (CallFlowandEMultipliers == true) QMult = ReadFlowandEMultipliers( "QMult",ISV, K, IAS, K, IAS);
				Tadv = 0;//initialize
				Tdisp = 0;//initialize
				
				//12_14_16 change
				if (Compartment[IAS][3] != 0) {
					Tdisp = EPrime * EMult; //IAS is boundary or dummy and WC has already been canceled
				}else{
					if (WCint == -999 ){
						System.out.println( " Error: WCint = -999 in Transport");
					}
					Tdisp = EPrime*EMult*WCint/WCk;	
				}

				
				if (Q == -999) {
					System.out.println("  Error: Q = -999");
					System.exit(1);
				}
				if (WCk == -999) {
					System.out.println("  Error: WCk = -999");
					System.exit(1);
				}
				if (Q < 0) Tadv = (float) (QMult*(Q*((float)1.0-AlphaToUse))/WCk);
				if (Q > 0) Tadv = (float) (QMult*Q*AlphaToUse/WCk);
				T = Tadv+Tdisp;
				if(VerboseLevel == 3) System.out.println("T = " +T);
				if(VerboseLevel == 3) System.out.println("Boundary condition type for adj compartment IAS = " +Compartment[IAS][3]);
				if (Compartment[IAS][3] == 0) { // Not a boundary or dummy compartment
					if(VerboseLevel == 3) System.out.println("Adj Compartment IAS not a boundary or dummy");
					if (SectnToEqn[ISV][IAS] != -999) { // ISV relevant in IAS
					// put transport terms in off-diagonal
					Row = SectnToEqn[ISV][K];
					Col = SectnToEqn[ISV][IAS];
					Asm[Row][Col] = T;
					if(VerboseLevel == 3) System.out.println("Asm["+Row+"]["+Col+"] = "+Asm[Row][Col]);
					} // end if (Compartment[IAS][3] == 0)
				} // end if
				
				if (Compartment[IAS][3] == 2) { // conventional (fixed) boundary condition
					if(VerboseLevel == 3) System.out.println("Adj Compartment IAS is a fixed boundary condition compartment");
					// Put transport terms in Bsm vector
					if (SectnToEqn[ISV][K] == -999) {
						System.out.println("  Error: SectnToEqn[ISV][K] = -999");
						System.exit(1);
					}
					
					JStart = 1; // initialize
					if (ISV == 1) {
						JStart = 1;
					}
					else {
						JStart += NT;
					} // end if
					for (JJ = JStart; JJ < MaxNSVxNT+1; JJ++) { // WALK BOUNDARY ARRAY TO FIND COMPARTMENT IAS for ISV
						if (Boundary[JJ][2] == IAS) {
							break; // Row JJ corresponds to boundary compartment IAS
						}
					} // end for JJ
					Bsm[SectnToEqn[ISV][K]] += T*Boundary[JJ][NColBoundary]; // units are gm/day
						//NOTE -- FOR CASE IS = 3, THE ZERO GRADIENT BOUNDARY CONDITION, THERE ARE NO OFF-DIAGONAL ELEMENTS
						//BECAUSE Sj = Si ALL "j'S" TERMS GET PUT INTO THE MAIN DIAGONAL AS i TERMS
	
				} // end if (Compartment[IAS][3] == 2)
				
				//MB code
				if(DoMassBalance == true){
					if(Compartment[IAS][3] != 0){
					  MBBoundaryRowCount += 1;
					  MBBoundary[MBBoundaryRowCount][1] = ISV;
					  MBBoundary[MBBoundaryRowCount][2] = K;
					  MBBoundary[MBBoundaryRowCount][3] = IAS;
					  if(Compartment[IAS][3] > 0 && Compartment[IAS][3] != 3){// dummy, conventional, or linear gradient boundary compartment
					  	  MBBoundary[MBBoundaryRowCount][4] = Tdisp * WCk + Tadv * WCk;// coefficient of c in off-diagonal, i.e. Cj in mb
					  }
					  if(Compartment[IAS][3] == 3){// zero gradient
					  	  MBBoundary[MBBoundaryRowCount][4] = 0;// coefficient of c in off-diagonal, i.e. Cj in mb
					  }
					}
				}
				// end off-diagonal

				// T WILL NOW BE TRANSPORT TERMS ASSOCIATED WITH Compartment K, I.E. THE MAIN DIAGONAL TERMS
				EMult = (float) 1.0;  //default
				QMult = (float) 1.0;  //default
				if (CallFlowandEMultipliers == true) EMult = ReadFlowandEMultipliers( "EMult",ISV, K, K, K, IAS);
				if (CallFlowandEMultipliers == true) QMult = ReadFlowandEMultipliers( "QMult",ISV, K, K, K, IAS);
	
				if (Compartment[IAS][3] == 3) { // zero gradient BC
					if (SectnToEqn[ISV][K] == -999) {
						System.out.println("  Error: SectnToEqn[IAS][K] = -999");
						System.exit(1);
					}
					Row = SectnToEqn[ISV][K];
					Col = SectnToEqn[ISV][K];
					Asm[Row][Col] += Q*QMult/WCk;
					if(VerboseLevel == 3) System.out.println("Asm["+Row+"]["+Col+"] = "+Asm[Row][Col]);
				}
				else { // all other COMPARTMENT types
					
					//12_14_16 change
				if (Compartment[IAS][3] == 0) {
					if (WCint == -999 ){
						System.out.println( " Error: WCint = -999 in Transport");
					}
					Tdisp = -EPrime*EMult*WCint/WCk;//IAS is interior
				}else{
					Tdisp = -EPrime*EMult; //IAS is not interior and old way is OK	
				}
					//Tdisp = -EPrime*EMult;
					

					if (Q < 0) Tadv = (float) QMult*Q*AlphaToUse/WCk;
					if (Q > 0) Tadv = (float) QMult*(Q*((float) 1.0-AlphaToUse))/WCk;
					T = Tadv+Tdisp;
					if(VerboseLevel == 3) System.out.println("T = " +T);
					if (SectnToEqn[ISV][K] == -999) {
						System.out.println("  Error: SectnToEqn = -999");
						System.exit(1);
					}
					Row = SectnToEqn[ISV][K];
					Col = SectnToEqn[ISV][K];
					Asm[Row][Col] += T;
					if(VerboseLevel == 3) System.out.println("Asm["+Row+"]["+Col+"] = "+Asm[Row][Col]);
					if(VerboseLevel == 3) Asm[SectnToEqn[ISV][K]][SectnToEqn[ISV][K]] += T;
 				} // end if (Compartment[IAS][3] == 3)
				//MB code
				if (DoMassBalance == true) {
					if(Compartment[IAS][3] > 0 && Compartment[IAS][3] != 3) {//dummy, conventional, or linear gradient boundary compartment
						MBBoundary[MBBoundaryRowCount][5] = Tdisp * WCk + Tadv * WCk;//coefficient of c in main-diagonal, i.e. Ci in mb	
					}
					if(Compartment[IAS][3] == 3) {//zero gradient boundary compartment
						MBBoundary[MBBoundaryRowCount][5] = Q*QMult;//coefficient of c in main-diagonal, i.e. Ci in mb	
					}
				}				
				// end main diagonal
	
				if (Compartment[IAS][3] == 4) { // LINEAR GRADIENT BOUNDARY -- BOTH MAIN DIAG AND OFF DIAG ELEMENTS NOW GET MODIFIED
					//NUMBER OF CompartmentS (INCLUDING BOUNDARY) ADJACENT TO Compartment K IS NAS
					//Get length (Avgl) between K and IAS
					 Avgl = GetLength(1, K, IAS);// use X,Y,Z coordinates
					 if(VerboseLevel == 3) System.out.println("Avgl = "+Avgl);
					//FIND NUMBER OF INTERIOR ADJ CompartmentS (NASINT) AND AVG LENGTH (AVGADJL)OVER THESE
					NASINT = 0;
					AvgAdjL = 0;
					JC = 2; // initialize to satisfy compiler
					AvgAdjL = 0; //initialize to satisfy compiler
					for (IK = 1; IK < NAS+1; IK++){
						if (IK == 1) {
							JC = 2;
						} else {
							JC +=  NColPerInface;
						} // end if
						IADJ = (int) Inface[K][JC]; // adjacent compartment number
						if (IADJ == -999) {
							System.out.println("  Error: IADJ = -999");
							System.exit(1);
						}
						if (Compartment[IADJ][3] ==0) { // interior compartment
							NASINT += 1;
							AvglNew = GetLength(1, K, IADJ);///only X,Y,Z coordinate method used
							if(VerboseLevel == 3) System.out.println("AvglNew = "+AvglNew);
							AvgAdjL += AvglNew;
						} //end if (Compartment[IADJ][3] == 0
					} // end for IK loop
					
					if (NASINT > NAS) {
						System.out.println("  Error: NASINT > NAS");
						System.exit(1);
					}
					if (NASINT == 0) {
						System.out.println("  Error: NASINT = 0");
						System.exit(1);
					}
					AvgAdjL = AvgAdjL/NASINT;
	
					// main diagonal terms
					if (SectnToEqn[ISV][K] == -999) {
						System.out.println("  Error: SectnToEqn[ISV][K] = -999");
						System.exit(1);
					}
					Row = SectnToEqn[ISV][K];
					Col = SectnToEqn[ISV][K];
					Asm[Row][Col] += T*(1.0+Avgl/AvgAdjL);
	
					// off-diagonal terms
					for (IK = 1; IK < NAS+1; IK++){
						if (IK == 1) {
							JC = 2;
						}
						else {
							JC += NColPerInface;
						} // end if
						IADJ = (int) Inface[K][JC]; // adjacent compartment number
						if (IADJ == -999) {
							System.out.println("  Error: IADJ = -999");
							System.exit(1);
						}
						if (Compartment[IADJ][3] == 0) { // interior compartment
							if (AvgAdjL*NASINT == 0) {
								System.out.println("  Error: AvgAdjL*NASINT  = 0");
								System.exit(1);
							}
							if (SectnToEqn[ISV][K] == -999) {
								System.out.println("  Error: SectnToEqn[ISV][K] == -999");
								System.exit(1);
							}
							Row = SectnToEqn[ISV][K];
							Col = SectnToEqn[ISV][IADJ];
							Asm[Row][Col]  -= T*Avgl/(AvgAdjL*NASINT); 
							//MB code
                                                        if(DoMassBalance == true){
                                                                MBBoundary[MBBoundaryRowCount][4] -= T * Avgl / (AvgAdjL * NASINT) * WCk;//coefficient of c in off-diagonal, i.e. Cj in mb
                                                        }
						} // end if (Compartment[IADJ][3] == 0
					} // end IK loop
					
				}// end if (Compartment[IAS][3] == 4)
			} // end J adj compartment loop	
			if (FlowTolerance != -999){
				float balance = Math.abs(100*(SumInflow - SumOutflow)/SumInflow);
				if (balance > FlowTolerance) {
					System.out.println("  Error: Flow balance > FlowTolerance in Transport()");
					System.exit(1); 
				}
			}
	
		} // end Compartment[K][3] == 0 && SVComp[ISV][K] == 1 if
	} //end for K loop
		
} // end Transport() method

 static float GetLength(int Option, int I, int J) {
	if(VerboseLevel > 1) System.out.println("I'm in GetLength() method ");
	if(VerboseLevel == 3) System.out.println("Option = "+Option+" Compartments "+I+ " and "+J);
	//following comment added by KWL 3/5/12
	//Note:  Option == 2 will never occur in this method because MixLengthOption = 2 is intercepted elsewhere.  Thus, the centroid-based data in "else if" below are never used.  This is true in VB version also.
	float Length = -999;
	if (Option == 1) {
			// COMPUTE LENGTH BETWEEN Compartment I AND J USING CENTROID COORDS
			// Compartment(*,4) IS X COORD; Compartment(*,5) IS  Y COORD; Compartment(*,6) IS Z COORD
			if(Compartment[I][4] == -999 || Compartment[J][4]== -999){
				System.out.println("Error:  Compartment[I][4] == -999 || Compartment[J][4}== -999 in GetLength method");
				System.exit(1);
			}
			if(Compartment[I][5] == -999 || Compartment[J][5]== -999){
				System.out.println("Error:  Compartment[I][5] == -999 || Compartment[J][5}== -999 in GetLength method");
				System.exit(1);
			}
			if(Compartment[I][6] == -999 || Compartment[J][6]== -999){
				System.out.println("Error:  Compartment[I][6] == -999 || Compartment[J][6}== -999 in GetLength method");
				System.exit(1);
			}			
			Length = (float) Math.pow((Compartment[I][4] - Compartment[J][4]),2);
			Length +=  (float) Math.pow((Compartment[I] [5] - Compartment[J][5]),2);
			Length +=  (float) Math.pow((Compartment[I][6] - Compartment[J][6]),2);
			Length = (float) Math.pow(Length,0.5);
	   if(VerboseLevel == 3) System.out.println("Length = " +Length);     	
	} else if(Option == 2) {
		System.out.println("  VB code not translated yet in GetLength() method !! ");
		System.exit(1); // error
	/*
	       Dim NumAdjSectn As Integer
		Dim AdjSectnNum As Integer
		Dim K As Integer
		Dim IC1 As Integer
		NumAdjSectn = Compartment(I, 2)
		For K = 1 To NumAdjSectn 'walk adjacent compartments
		    If K = 1 Then
			IC1 = 2
		    Else
			IC1 = IC1 + NColPerInface
		    End If
		    AdjSectnNum = Inface(I, IC1)
		    If AdjSectnNum = I Then Stop 'data error
		    Length = Inface(I, IC1 + 2) 'distance from centroid to interface with AdjSectnNum
		Next K
		'get length from centroid J to interface with I
		NumAdjSectn = Compartment(J, 2)
		For K = 1 To NumAdjSectn 'walk adjacent compartments
		    If K = 1 Then
			IC1 = 2
		    Else
			IC1 = IC1 + NColPerInface
		    End If
		    AdjSectnNum = Inface(J, IC1)
		    If AdjSectnNum = J Then Stop 'data error
		    Length = Length + Inface(J, IC1 + 2) 'distance from centroid I to centroid J
		Next K
	*/
	} else {
		System.out.println("  Error: in GetLength, Option neither 1 nor 2");
		System.exit(1); // error
	} // end if
	if (Length == -999) {
		System.out.println("Error!  Length = -999 in GetLength");
		System.exit(1);
	}
	return Length;

} // end GetLength() method

	static float ReadFlowandEMultipliers(String Which,int ISV, int CompRow, int CompCol , int IntI , int IntJ){
		float Multiplier = (float)1.0;
		if(VerboseLevel > 1) System.out.println("I'm in ReadFlowandEMultipliers() ");
		//field 1 is state variable, field 2 is compartment I, field 3 is compartment J, field 4 is I,J E multiplier, field 5 is I,J Q multiplier
		int SV2 = -999;
		try{
			String strFile = "FlowandEMultipliers.csv";
			BufferedReader br = new BufferedReader( new FileReader(strFile));//create BufferedReader to read csv file
			String strLine = "";
			StringTokenizer st = null;
			strLine = br.readLine();//SV,Equation Compartment ,Compartment Containing Multiplier,Interface I,Interface J,Eij multiplier, Qij multiplier
			if(VerboseLevel == 3) System.out.println(strLine);
		
			do {
				strLine = br.readLine();//data record
				if(VerboseLevel == 3) System.out.println(strLine);
				st = new StringTokenizer(strLine, ",");
				
				SV2 = Integer.parseInt(st.nextToken());
				if (SV2 == -999) break;// end of file
				//CallFlowandEMultipliers = true;  deleted line by KWL 3/5/12
				int I = Integer.parseInt(st.nextToken());
				int J = Integer.parseInt(st.nextToken());
				int K = Integer.parseInt(st.nextToken());
				int L = Integer.parseInt(st.nextToken());
				float EMult2 = Float.parseFloat(st.nextToken());
				float QMult2 = Float.parseFloat(st.nextToken());
				if (I == CompRow && J == CompCol){
					if ((K == IntI && L == IntJ)||(K == IntJ && L == IntI)){
                                            if(Which.equals("EMult")){
						Multiplier = EMult2;
					    }else{
						Multiplier = QMult2;
                                        }
						break; // exit do
					}
				}
		
				} while (SV2 != -999  );
			br.close(); // close input stream
		} // end try
		catch(Exception e)
		{
			System.out.println("Exception while reading FlowandEMultipliers.csv file: " + e); 
		} // end catch
	     return Multiplier;
	} // end ReadFlowandEMultipliers method
	
	static void ForceFn(float Time) {
		if(VerboseLevel> 1) System.out.println("I'm in ForceFn()");
		if(VerboseLevel == 3) System.out.println("Time = "+Time);
		float TimeStart, TimeStop, Load;
		int StartRow = 0;
		int dummycounter = 0;
		int Comp,  Row, Counter,NumDataRecords, I, J,SV;
		String dummy;
		
		if (EquationSolver == false) {
			for ( I = 1; I < NEQN + 1; I++){
				b[I] = (float) 0.0; // initialize
				//MB code
				if(DoMassBalance == true){
					MBLoad[I]=0;//initialize	
				}
			} //end for I loop
		
			try
			{
				String strFile = "Loads.csv";
				BufferedReader br = new BufferedReader( new FileReader(strFile));//create BufferedReader to read csv file
				String strLine = "";
				StringTokenizer st = null;
				do {// time block loop
					strLine = br.readLine();//header record: TimeStart, TimeEnd
					if(VerboseLevel == 3) System.out.println("strLine = " +strLine);
					
					if (strLine.equals("-999")) {
						if(VerboseLevel == 3) System.out.println("exiting time block loop");
						break; // end of file -- exit do
					}
					strLine = br.readLine();
					st = new StringTokenizer(strLine, ",");
					if(VerboseLevel == 3) System.out.println("strLine = "+strLine);
					TimeStart = Float.parseFloat(st.nextToken());
					if(VerboseLevel == 3) System.out.println("TimeStart = "+ TimeStart);
					TimeStop = Float.parseFloat(st.nextToken());
					if(VerboseLevel == 3) System.out.println("TimeStop = "+ TimeStop);
					if (Time >= TimeStart && Time <= TimeStop) {
						if(VerboseLevel == 3) System.out.println("top of time block if statement");
						strLine = br.readLine();//Compartment,SV, SV load, SV, SV load, ...
						if(VerboseLevel == 3) System.out.println("strLine = "+strLine);
						do{// loop over compartments within time block
							strLine = br.readLine();// data record
							st = new StringTokenizer(strLine, ",");
							Comp = Integer.parseInt(st.nextToken());
							if(VerboseLevel == 3) System.out.println("Comp = "+Comp);
							if (Comp == -999) {
								if(VerboseLevel == 3) System.out.println("exiting do/while compartment loop");
								break; // end of file -- exit do
							} // end if 
							do{ // loop over SVs and SV loads within a compartment
								SV = Integer.parseInt(st.nextToken());
								//System.out.println("SV = "+SV);
								if (SV == -999) {
									if(VerboseLevel == 3) System.out.println("exiting do/while SV within compartment loop");
									break; // exit do
								} // end if 		
								Load = Float.parseFloat(st.nextToken());
								if(VerboseLevel == 3) System.out.println("load = " +Load);
								if (Comp == 1){
									StartRow=0;
								}else{
									StartRow += NEQNsm[SV-1];
								}
								Row = (int) StartRow + SectnToEqn[SV][Comp];
								b[Row] = Load/Compartment[Comp][9];//divide by water content
							        //MB code
							        if(DoMassBalance == true){
							        	MBLoad[Row] = Load;	
							        }
								if (Linear == false && Explicit == false){
									b[Row] = b[Row] * NewtonSF;
								} // end if (Linear == false && Explicit == false)
							} while (dummycounter == 0);// end loop over SVs and SV loads within compartment
						}while(dummycounter == 0);// end loop over compartments within time block
						if (Comp == -999){
							if(VerboseLevel == 3) System.out.println("exiting time block loop");
							break;
						}
					} else {
						do {// jumping over not-applicable time blocks
							if(VerboseLevel == 3) System.out.println("reading...");
							strLine = br.readLine();
							if(VerboseLevel == 3) System.out.println("strLine = " +strLine);
							if (strLine.equals("-999")) {
								if(VerboseLevel == 3) System.out.println("ending non-applicable time block loop");
								break; // end of file -- exit do
							}
						}while(dummycounter == 0);//end jumping over non-applicable time blocks loop
						if (Dynamic == false) break; // exit do -- read only one time block for steady-state
					} // end if (Time >= TimeStart && Time <= TimeStop)
				} while(dummycounter == 0); // end next timestart, timestop block do
				if(VerboseLevel == 3) System.out.println("exiting reading forcing fn");
				br.close(); // close input stream
			}// end try
			catch(Exception e) {
				System.out.println("Exception while reading Loads.csv file: " + e); 
			}// end catch
			// Now add boundary condition terms
			for (I = 1;I < NEQN+1;I++){
				if(Dynamic == false && Explicit == false){
					BBC[I]=BBC[I]*NewtonSF;
				}
				b[I]=b[I]+BBC[I];// on LHS
				//start KWL 3/5/12 change
				//If Dynamic = False Then B(I) = -B(I) 'ON RHS
				if(Dynamic == false) b[I]=-b[I];
				//end KWL 3/5/12 change
			}//end I loop
		}else{ // EquationSolver = true
		    if(ShellBEqnSolver== true) {
			WriteTime("BEqnSolver",Time);
			ShellFunction("BEqnSolver.exe");
		    }
		    try {
		       String strFile = externalBfile;
		       BufferedReader br = new BufferedReader( new FileReader(strFile));//create BufferedReader to read csv file
		       String strLine = "";
		       StringTokenizer st = null;
		        
			 for (int i=1; i <NEQN+1;i++) {
				 strLine = br.readLine();
				 b[i] = Double.parseDouble(strLine);
				 if(Linear == false && Explicit == false) b[i]=b[i]*NewtonSF;
				 if(VerboseLevel == 3) System.out.println("b["+i+"] = " +b[i]);
			 } // end for i
			 br.close(); // close input stream
		    }// end try
		    catch(Exception e)
			 {
				System.out.println("Exception while reading "+e);
				System.exit(1);
			 }// end catch

		}// end if (EquationSolver == false)
	} // end ForceFn() method
//start KWL 3/5/12 changes to LinearSolver (whole thing is changed)
static void LinearSolver( int n,  double A[][],  double b[])
	{
	/*
	This algorithm was found on the internet and is under the GNU Public License.
	It was developed by Jon Squire, an Adjunct Faculty in Computer Science and Electrical Engineering 
	at the University of Maryland Baltimore County (UMBC).
	It uses Gauss-Jordan elimination to reduce the "B" matrix, which has the RHS as the last column, to the 
	identity matrix, so that the solution can be directly read from the last (RHS) column.  It uses pivoting
	on the largest diagonal element to avoid division by zero and to improve numerical accuracy
	A discussion and example can be seen at: http://www.cs.umbc.edu/~squire/cs455_lect.html#L3
	It has been modified by K Little for the GEM project to "ignore the 0 element" in arrays, to 
	be consistent with the rest of the GEM array treatment.  It is also modified to use the A and b arrays without cramming b into A as the last column.
	*/

	if(VerboseLevel> 1) System.out.println("I'm in LinearSolver ");
	if(VerboseLevel == 3) System.out.println("length of A array = "+A.length);
	if(VerboseLevel == 3) System.out.println("length of b array = "+b.length);
	
	int hold, I_pivot,I,J;
	double pivot;
	double abs_pivot;
	int row[] = new int[n+1];
	double X[] = new double[n+1];
	
	if(VerboseLevel == 3){
		for(I=1;I<n+1;I++){
			 System.out.println("b["+I+"] = "+b[I]);
			for (J=1;J<n+1;J++){
				System.out.println("A["+I+"]["+J+"] = "+A[I][J]);
			}
		}
	}
	
	//set up row interchange vectors
	for(int k=1;k<n+1;k++){
		row[k]=k;
	}
	//main reduction loop
	if(VerboseLevel == 3) System.out.println("main reduction loop ...");
	for (int k=1;k<n+1;k++){
		//find largest pivot element
		if(VerboseLevel == 3) System.out.println("k = "+k+"and finding largest pivot element");
		pivot = A[row[k]][k];
		abs_pivot=Math.abs(pivot);
		I_pivot = k;
		for(int i=k;i<n+1;i++){
			if(Math.abs(A[row[i]][k]) > abs_pivot){
				I_pivot = i;
				pivot = A[row[i]][k];
				abs_pivot = Math.abs(pivot);
			}
		}
		//have pivot, interchange row indices
		if(VerboseLevel == 3) System.out.println("pivot,interchange row indices ...");
		hold = row[k];
		row[k]=row[I_pivot];
		row[I_pivot]=hold;
		//check for near-singular
		if(VerboseLevel == 3) System.out.println("check for near-singularity...");
		if(abs_pivot < 1.0E-10){
			for(int j=k+1;j<n+2;j++){
				if (j < n+1){
					A[row[k]][j] = 0.0;
				}else{
					b[row[k]]=0;
				}
			}
			System.out.println("redundant row (singular) "+row[k]);
		}else{
			//reduce about pivot
			if(VerboseLevel == 3) System.out.println("reduce about pivot ...");
			for(int j=k+1;j<n+2;j++){
				if(VerboseLevel == 3) System.out.println("j = "+j);
				if(VerboseLevel == 3) System.out.println("row[k] = "+row[k]);
				//if(VerboseLevel == 3) System.out.println("A["+row[k]+"]["+j+"] = "+A[row[k]][j]);
				if(VerboseLevel == 3) System.out.println("A["+row[k]+"]["+k+"] = "+A[row[k]][k]);
				if(j<n+1){
					A[row[k]][j] = A[row[k]][j]/A[row[k]][k];
				}else{
					b[row[k]]=b[row[k]]/A[row[k]][k];
				}
			}
			//inner reduction loop
			if(VerboseLevel == 3) System.out.println("inner reduction loop");
			for(int i=1;i<n+1;i++){
				if(i!=k){
					for(int j=k+1;j<n+2;j++){
						if(j<n+1){
							A[row[i]][j]=A[row[i]][j]-A[row[i]][k]*A[row[k]][j];
						}else{
							b[row[i]]=b[row[i]]-A[row[i]][k]*b[row[k]];
						}
					}
				}
			}//finished inner reduction

		} //end main reduction loop
	}//end k loop		
	//build X for return, unscrambling rows
	for(int i=1;i<n+1;i++){
		X[i]=b[row[i]];
	}
	for(int i=1;i<n+1;i++){
		b[i]=X[i];
	}	

	}// end LinearSolver method
//end KWL 3/5/12 changes to LinearSolver
static void EulerOption(int TimeStep, float CumTime, boolean OutputTime, float TempTime) {
	
		int I,J;
		String EulerStabilityfname = "EulerStabilityInfo.dng";
		if(TimeStep == 1){
			if(VerboseLevel > 1) System.out.println("Implementing Euler option ...");
		}
		if(VerboseLevel> 1) System.out.println("I'm in EulerOption with TimeStep = " + TimeStep+" and CumTime = "+CumTime);
		if(VerboseLevel == 3) System.out.println("Linear = " +Linear);
		if(VerboseLevel == 3) System.out.println("Explicit = " +Explicit);
		if (Linear == true || Explicit == true){

			//c^t+1 = ((RV)^-1)^t+1 [(RVc)^t + delT(Ac + b)^t]
			if (TimeStep == 1) {
				AC = new double[NEQN+1];
				RVt= new double[NEQN+1];
				RVtp1= new double[NEQN+1];
			}
			if(VerboseLevel == 3) System.out.println("DynParm = " +DynParm);
			System.out.println("here and DynParm = "+ DynParm);
			if (DynParm == true){
				ParmUpdate(CumTime);
			} 	
			
			for (I = 1; I < NEQN+1; I++){
				if (EquationSolver == false) {
					RVt[I] = R[I]*V[I]; // R*V at t
					if(VerboseLevel == 3) System.out.println(" at I = "+I);
					if(VerboseLevel == 3) System.out.println("R = " +R[I]);
					if(VerboseLevel == 3) System.out.println("V = " +V[I]);
				} else{
					RVt[I] = RVEqnSolver[I];
				}
			}// end I loop
			
			if (EquationSolver == false){
				//load transport terms into array A 
				LoadATrans();
								
				// load linear source/sink terms into array A 
				LoadALinearSrcSnk();
								
				if (TimeStep == 1){// create diagnostics file
					try {
						FileWriter fw = new FileWriter(EulerStabilityfname);
						BufferedWriter out = new BufferedWriter(fw);
						//out.write("Euler Option stability information");
						out.close();
					} catch (IOException e){
					 	 e.printStackTrace();
					}
				}
				if (CallWiggle == true){
					Wiggle(CumTime);
				}
         
			}// end if (EquationSolver == false)
			ForceFn(CumTime);
			if (OutputTime == true){
				TempTime = CumTime;
				Aout(TempTime);// KWL 3/14 change
				RVout(TempTime);// KWL 3/14 change
				Bout(TempTime);// KWL 3/14 change
			}

			

                        for (I = 1; I < NEQN+1;I++){
                        	AC[I] = 0; // initialize
                        	for (J = 1; J < NEQN+1; J++){
                        		AC[I] += A[I][J]* c[J];
                        	}
                        } 
                        	
                        if (Linear == false){
                        	FofC(c,TempTime);// evaluate at timestep 't' for Explicit == true
                        }
                        
                        if (DynParm == false){
                        	for (I = 1; I < NEQN+1; I++){
                        		if (EquationSolver == false){
                        			RVtp1[I] = RVt[I]; // R*V at t+1
                        		} else{
                        			RVtp1[I] = RVEqnSolver[I];
                        		}
                        	}
                        	
                        } else{ // DynParm = true
                        	ParmUpdate(CumTime+DelT);// need volume, retardation at t+1
                        	for (I = 1; I < NEQN+1; I++){
                        		if (EquationSolver == false){
                        			RVtp1[I] = R[I]*V[I]; // R*V at t+1
                        		} else{
                        			RVtp1[I] = RVEqnSolver[I];
                        		}
                        	}
                        }// end else
                        for (I = 1; I<NEQN+1;I++){
                        	if (EquationSolver == false){
                        		if (Linear == true){
                        			if(VerboseLevel == 3) System.out.println("RVTp1["+I+"]="+RVtp1[I]);
                        			if(VerboseLevel == 3) System.out.println("RVT["+I+"]="+RVt[I]);
                        			if(VerboseLevel == 3) System.out.println("b["+I+"]="+b[I]);
                        			if(VerboseLevel == 3) System.out.println("Delt="+DelT);
                        			if(VerboseLevel == 3) System.out.println("AC["+I+"]="+AC[I]);
                        			if(VerboseLevel == 3) System.out.println("before update c["+I+"] = "+c[I]);  
                        			if (DynParm == true || TimeStep == 1){
                                         		StableCriterion = RVt[I]+DelT*A[I][I];
                        				if (StableCriterion < 0){
								DTmax = -RVt[I]/A[I][I];
								try {
									BufferedWriter out = new BufferedWriter(new FileWriter(EulerStabilityfname, true));
                        						//out.newLine();
                        						out.write("\r\n");// line feed
                        						out.write("stability criterion violated at time "+CumTime+" for equation "+I);
                        						out.write(" max time step is "+DTmax);	
                        						out.close();
								} catch (IOException e){
                        						e.printStackTrace();
                        					}
							}
                        			}
                        			c[I] = (1/RVtp1[I])*(RVt[I]*c[I]+DelT*(AC[I]+b[I]));
                        		} else{ // nonlinear
                        			if (DynParm == true || TimeStep == 1){
                        				StableCriterion = RVt[I]+DelT*(A[I][I]+f[I]);
                        				System.out.println("StableCriterion = "+StableCriterion);
                        				if (StableCriterion < 0){
								DTmax = -RVt[I]/(A[I][I]+f[I]);
								try {
									BufferedWriter out = new BufferedWriter(new FileWriter(EulerStabilityfname, true));
                        						//out.newLine();
                        						out.write("\r\n");// line feed
                        						out.write("stability criterion violated at time "+CumTime+" for equation "+I);
                        						out.write(" max time step is "+DTmax);	
                        						out.close();
								} catch (IOException e){
                        						e.printStackTrace();
                        					}
							}
						}
                         			c[I] = (1/RVtp1[I])*(RVt[I]*c[I]+DelT*(AC[I]+b[I]+f[I]));
                        		}
                        		} else { // EquationSolver = true
                        			System.out.println("code needed here in EulerOption");
                        			/*
					       If Linear = True Then
					       	c(I) = (1 / RVEqnSolver(I)) * (RVEqnSolver(I) * c(I) + DelT * (AC(I) + B(I)))
					       	Else ' nonlinear
					       		c(I) = (1 / RVEqnSolver(I)) * (RVEqnSolver(I) * c(I) + DelT * (AC(I) + B(I) + f(I)))
						End If
                        		*/
                        	}
                        	if(VerboseLevel == 3) System.out.println("** c["+I+"] = "+c[I]);                     	
                        } // end I loop
		} else { 
			//start KWL 3/14 changes
			Newton(cupdate,c,CumTime);
		
			for(I = 1;I < NEQN+1;I++){
				c[I] = cupdate[I];
			}
			//end KWL 3/14 changes
		}// end if (Linear == true || Explicit == true)
	}// end EulerOption()
	
	static void MacCormack(int TimeStep, float CumTime, boolean OutputTime, float TempTime){
		int I,J;
		if(EquationSolver == true || Linear == false) {
			System.out.println("MacCormack method not applicable -- choose another dynamic method");
			System.exit(1); 
			}
		if(Rnonlinear == false) {
			//c^t+1 = [(RV)^-1)^t+1 (RVc)^t] + (delT/2){[A(alpha = 0)c + f(c) + b]^t + [A(alpha = 1)chat + f(chat) + b]^t+1}
			if(TimeStep == 1){
				if(VerboseLevel > 1) System.out.println("Implementing MacCormack option ...");
				Slope1 = new double[NEQN+1];
				Slope2 = new double[NEQN+1];
				chat = new double[NEQN+1];
				RVt = new double[NEQN+1];
				RVtp1 = new double[NEQN+1];
			}
			if(DynParm == true){
				ParmUpdate(CumTime);// time t
			}
			for (I = 1; I < NEQN+1; I++){
				RVt[I] = R[I]*V[I]; // R*V at t
				}
				// calculate Slope1 and Chat at time 't'
				Alpha = 0;//forward difference
				//load transport terms into array A 
				LoadATrans();
				// load linear source/sink terms into array A 
				LoadALinearSrcSnk();
				if (CallWiggle == true){
					Wiggle(CumTime);
				}
				ForceFn(CumTime);// time t
				if (OutputTime == true){
					TempTime = CumTime;
					Aout(TempTime);// KWL 3/14 change
					RVout(TempTime);// KWL 3/14 change
					Bout(TempTime);// KWL 3/14 change
			}
			if (Linear == false){
				FofC(c,TempTime);// evaluate at timestep 't' for Explicit == true
			}
			for (I = 1; I < NEQN+1;I++){
					double SumAC = 0;
					for (J = 1; J < NEQN+1; J++){
						SumAC += A[I][J]* c[J];
					}
					Slope1[I] = SumAC+b[I];
					if(Linear == false){
						Slope1[I] += f[I];
					}
			} 
			if(DynParm == false){
				for(I=1;I<NEQN+1;I++){
					RVtp1[I]=RVt[I];//R*V at t+1
				}
			} else{
				ParmUpdate(CumTime+DelT);//need volume, retardation at t+1
				for(I=1;I<NEQN+1;I++){
					RVtp1[I]=R[I]*V[I];//R*V at t+1
				}
			}
			for(I=1;I<NEQN+1;I++){
				chat[I]=(1/(RVtp1[I]))*(RVt[I]*c[I]+DelT*Slope1[I]);
			}
			
			//calculate Slope2 using Chat at time 't+1'
			if(DynParm == true){
				ParmUpdate(CumTime+DelT);// already have RVt and RVtP1 from above
			}
			Alpha = 1;//backward difference
			LoadATrans();
			// load linear source/sink terms into array A 
			LoadALinearSrcSnk();
			//get loads at t+1
			ForceFn(CumTime+DelT);
			if(Linear == false){
				TempTime = CumTime + DelT;
				FofC(chat,TempTime);
			}
			for (I = 1; I < NEQN+1;I++){
					double SumAC = 0;
					for (J = 1; J < NEQN+1; J++){
						SumAC += A[I][J]* chat[J];
					}
					Slope2[I] = SumAC+b[I];
					if(Linear == false){
						Slope2[I] += f[I];
					}
			} 
			for(I=1;I<NEQN+1;I++){
				c[I]=(1/(RVtp1[I]))*(RVt[I]*c[I]+DelT*(Slope1[I]+Slope2[I])/2);
			}
		}else{
			System.out.println("Nonlinearnot supported by MacCormack method");
			System.exit(1);
		}//end if(Rnonlinear == false)
	} // end MacCormack method

	 static void BackTimeOption(int TimeStep, float CumTime, boolean OutputTime, float TempTime) {
	int I,J;
	if(VerboseLevel > 1) System.out.println("I'm in BackTimeOption with TimeStep = " + TimeStep+" and CumTime = "+CumTime);
	if(TimeStep == 1){
		if(VerboseLevel> 1) System.out.println("Implementing Back Time option ...");
	}
	
	if (Linear == true){
		//((RV)^t+1-(delT)A^t+1)c^t+1 = (RVc)^t + (delT)b^t+1
		if(TimeStep == 1){
			Aimplicit = new double[NEQN+1][NEQN+1];	
			Bimplicit = new double[NEQN+1];
			//Asolve = new double[NEQN+1][NEQN+2];//will hold RHS as last column for LinearSolver() method -- line changed by KWL 3/5/12
			// initialize
			for (I = 1; I < NEQN+1; I++){
				Bimplicit[I] = (float) 0.0;
				for (J = 1; J < NEQN+1; J++){
					Aimplicit[I][J] = (float) 0.0;
				} // end for J
			}// end for 	I
		}
		// SET UP MASS BALANCE EQUATIONS IN COEFFICIENT MATRIX A at time step 't+1'
		
		if(DynParm == true){
			ParmUpdate(CumTime+DelT);
		}
		if(EquationSolver == false){
			//load transport terms into array A
			LoadATrans();//BBC loaded with boundary conditions only and on LHS
			LoadALinearSrcSnk();
			if (CallWiggle == true){
				Wiggle(CumTime);
			}
	
		    if (OutputTime == true){
				TempTime = CumTime;
				Aout(TempTime);// KWL 3/14 change
				RVout(TempTime);//KWL 3/14 change
		    }
		}// end if(EquationSolver == false)
		for (I = 1; I < NEQN+1; I++){
				for (J = 1; J < NEQN+1; J++){
					if (I == J){
						if(EquationSolver == false){
							Aimplicit[I][J]=R[I]*V[I]-DelT*A[I][J];
						}else{
							//System.out.println("VB code needs translating here in BackTimeOption");
							//Aimplicit(I, J) = RVEqnSolver(I) - DelT * A(I, J)
						}
					}else{
						Aimplicit[I][J]= - DelT*A[I][J];
					}// end if I == J
				} // end for J
			}// end for I
	
		//get loads at TimeStep = 't+1'
		ForceFn(CumTime+DelT);//B returned with BC's and loads on correct side
		if (OutputTime == true){
			TempTime = CumTime+DelT;//KWL 3/15 change
		     Bout(TempTime);// KWL 3/14 change
	      }
	
	
		if(DynParm == true){
			ParmUpdate(CumTime);//get R and V at time = t
		}
		for (I = 1; I < NEQN+1; I++){
			if(EquationSolver == false){
				Bimplicit[I]= R[I]* V[I]* c[I]+ DelT * b[I];
			}else{
				System.out.println("VB code needs translating here in BackTimeOption");
				//Bimplicit(I) = RVEqnSolver(I) * c(I) + DelT * b(I)
			}
		}
		//start KWL 3/5/12 changes
		// load Asolve matrix for solve method putting b vector back to RHS
		//for(I=1;I<NEQN+1;I++){
			//for(J=1;J<NEQN+2;J++){
				//if(J<NEQN+1){
					//Asolve[I][J]=Aimplicit[I][J];
				//}else{
					//Asolve[I][J]=Bimplicit[I];
				//}
			//}
		//}
		//LinearSolver(NEQN,Asolve,c);
		LinearSolver(NEQN,Aimplicit,Bimplicit);
		for (I = 1; I < NEQN+1; I++){
			c[I] = Bimplicit[I];
		}
		//end KWL 3/5/12 changes
		if(VerboseLevel == 3){
			for (I = 1; I < NEQN+1; I++){
				System.out.println("c["+I+"]= "+c[I]);//updated for 't+1'
			}
		}
	}else{//nonlinear ...
		if (EquationSolver == false){
			if (CallWiggle == true){
				Wiggle(CumTime);
			}	
		}
		
		Newton(cupdate,c,CumTime);
		
		for(I = 1;I < NEQN+1;I++){
			c[I] = cupdate[I];
		}

		if(VerboseLevel == 3){
			for(I = 1;I < NEQN+1;I++){
				System.out.println("after call to Newton, cupdate["+I+"]="+cupdate[I]);
				System.out.println("after call to Newton, c["+I+"]="+c[I]);
			}
		}
		TempTime = CumTime+DelT;
	}// end if (Linear == true)
}//end BackTimeOption


static void CenterTimeOption(int TimeStep, float CumTime, boolean OutputTime, float TempTime) {
		int I,J;
		double SumAC;
		if(VerboseLevel> 1) System.out.println("I'm in CenterTimeOption with TimeStep = " + TimeStep+" and CumTime = "+CumTime);
		if(TimeStep == 1){
			if(VerboseLevel> 1) System.out.println("Implementing Center Time option ...");
		}
		
		if(Linear == true){
		     //((RV)^t+1-(delT/2)A^t+1)c^t+1 = (RVc)^t + (delT/2)(b^t+1 + (Ac)^t+1 + b^t+1)
		     if(TimeStep == 1){
		     	     Aimplicit = new double[NEQN+1][NEQN+1];	
		     	     Bimplicit = new double[NEQN+1];
		     	     //Asolve = new double[NEQN+1][NEQN+2];//will hold RHS as last column for LinearSolver() method -- KWL 3/5/12 change
		     	     // initialize
		     	     for (I = 1; I < NEQN+1; I++){
		     	     	     Bimplicit[I] = (float) 0.0;
		     	     	     for (J = 1; J < NEQN+1; J++){
		     	     	     	     Aimplicit[I][J] = (float) 0.0;
		     	     	     } // end for J
		     	     }// end for I
		     }//end if(TimeStep == 1)
		     if(EquationSolver == false){
			 //load transport terms into array A
			 LoadATrans();//BBC loaded with boundary conditions only and on LHS
			 LoadALinearSrcSnk();
			 if (CallWiggle == true){
				 Wiggle(CumTime);
			 }
			 if (OutputTime == true){
				 TempTime = CumTime;
				 Aout(TempTime);// KWL 3/14 change
				 RVout(TempTime);//KWL 3/14 change
			 }	
		     }// end if(EquationSolver == false)
		     for (I = 1; I < NEQN+1; I++){
			for (J = 1; J < NEQN+1; J++){
				if (I == J){
					if(EquationSolver == false){
						Aimplicit[I][J]=R[I]*V[I]-(DelT/2)*A[I][J];
					}else{
						System.out.println("VB code needs translating here in BackTimeOption");
						//Aimplicit(I, J) = RVEqnSolver(I) - (DelT / 2) * A(I, J)  'unitless
					}
				}else{
					Aimplicit[I][J]= - (DelT/2)*A[I][J];
				}// end if I == J
				} // end for J
		     }// end for I   
            	 
		     //get loads at TimeStep = 't+1'
		     ForceFn(CumTime+DelT);//B returned with BC's and loads on correct side
		
		     if (OutputTime == true){
		     	TempTime = CumTime+DelT;// KWL 3/14 change
			Bout(TempTime);// KWL 3/14 change
		     }

		     for (I = 1; I < NEQN+1; I++){
			Bimplicit[I]= (DelT/2) * b[I];
		     }
		     
		     if(EquationSolver == false) {
			if(DynParm == true){
				ParmUpdate(CumTime);
				//load transport terms into array A
				LoadATrans();//BBC loaded with boundary conditions only and on LHS
				LoadALinearSrcSnk();   	
			}
			  
		     }//end if(EquationSolver == false)
		     
		     //get loads at TimeStep = 't'
		     ForceFn(CumTime);//B returned with BC's and loads on correct side
		     
                    for(I = 1;I<NEQN+1;I++){
                    	    SumAC = 0;
                    	    for(J=1;J<NEQN+1;J++){
                    	    	    SumAC += A[I][J]*c[J];
                    	    }//end J loop
                    	    if (EquationSolver == false) {
                    	    	    Bimplicit[I] += R[I]*V[I]*c[I]+(DelT/2)*(SumAC+b[I]);
                    	    }else{
                    	    	   System.out.println("VB code needs translating here in BackTimeOption"); 
                    	    	   // Bimplicit(I) = Bimplicit(I) + RVEqnSolver(I) * c(I) + (DelT / 2) * (SumAC + B(I))
                    	    }
                    	    
                    }//end I loop
                    
                    //start KWL 3/5/12 changes
                    // load Asolve matrix for solve method putting b vector back to RHS
                    //for(I=1;I<NEQN+1;I++){
                    	    //for(J=1;J<NEQN+2;J++){
                    	    	   // if(J<NEQN+1){
                    	    	    	    //Asolve[I][J]=Aimplicit[I][J];
                    	    	    //}else{
                    	    	    	    //Asolve[I][J]=Bimplicit[I];
                    	    	    //}
                    	    //}
                    //}

                    //if(VerboseLevel == 3) System.out.println("before LinearSolver, Asolve(NEQN,NEQN)= "+Asolve[NEQN][NEQN]);
                    //LinearSolver(NEQN,Asolve,c);
                    LinearSolver(NEQN,Aimplicit,Bimplicit);
                    for(I=1;I<NEQN+1;I++){
                    	    c[I] = Bimplicit[I];
                    }
                    
                    //end KWL 3/5/12 changes
                    if(VerboseLevel == 3) {
                    	    for (I = 1; I < NEQN+1; I++){
                    	    	    System.out.println("c["+I+"]= "+c[I]);//updated for 't+1'
                    	    }
                    }
                    //start KWL 3/9/12 changes
		}else{//nonlinear
			if (EquationSolver == false){
				if (CallWiggle == true){
					Wiggle(CumTime);
				}	
			}	
					
			Newton(cupdate,c,CumTime);
			
			for(I = 1;I < NEQN+1;I++){
				c[I] = cupdate[I];
			}
	
			if(VerboseLevel == 3){
				for(I = 1;I < NEQN+1;I++){
					System.out.println("after call to Newton, cupdate["+I+"]="+cupdate[I]);
					System.out.println("after call to Newton, c["+I+"]="+c[I]);
				}
			}
			
			TempTime = CumTime+DelT;
			//end KWL 3/9/12 changes
		}//end if(Linear == true)
		
}//end CenterTimeOption

 static void Wiggle(float Time) {
	if(VerboseLevel > 1) System.out.println("I'm in Wiggle()");
	//DETERMINE PECLET NUMBER FOR SOLUTION POSITIVITY/NO WIGGLE FOR CENTRAL DIFFERENCES AND E, Q > 0
	int I,J,Ipass,ICInface,ICFlow,ICEcoef,K,L,NAdj,RowCount,ISV;
	float Q,E,V,WC,Length,Test,AlphaToUse,Area;
	boolean Found;
	String Pos;
	for(ISV = 1;ISV<MaxNSV+1;ISV++){
		RowCount = (ISV-1)*NT;
		for(I=1;I<NT+1;I++){
			RowCount +=1;
			if(Compartment[I][3] == 0) {//I is interior compartment
				if(VerboseLevel == 3) System.out.println("Interior compartment = "+I);
				NAdj = (int) Compartment[I][2];//number of compartments adjacent to I
				if(VerboseLevel == 3) System.out.println("number of adj compartments = "+NAdj);
				ICInface = 2;//initialize
				ICFlow = 2;//initialize
				ICEcoef = 3;//initialize
				for(J =1;J<NAdj+1;J++){
					if(J > 1){
						ICInface += NColPerInface;
						ICFlow += GEM.NColPerFlow;
						ICEcoef += GEM.NColPerEcoef;
					}
					K = (int)Inface[I][ICInface];// ADJ Compartment NUMBER
					if(VerboseLevel == 3) System.out.println("adj compartment K = "+K);
					Q = Flow[I][ICFlow + 1];// cu m/day
					E = ECoef[RowCount][ICEcoef + 1];// sq m/DAY
					WC = Compartment[I][9];// water content
					Alpha = Inface[I][ICInface + 4];
					if(VerboseLevel == 3) System.out.println("Q="+Q+" E= "+E+" WC = "+WC);
					if(Q<0 && Alpha == 0.5){//wiggle not an issue for positive flow
						Length = GetLength(1, I, K);//X,Y,Z coordinates used	
						if(VerboseLevel == 3) System.out.println("K = "+K);
						Area = Inface[I][ICInface + 1];
						V = Math.abs(Q) / Area;// VELOCITY IN m/DAY
						//find alpha
						Found = false;
						AlphaToUse = Alpha;//initialize
						for(L=1;L<NAdj+1;L++){	
							if(VerboseLevel == 3) System.out.println("L = "+L);
							if(Alpha == 0.5){
								if(VerboseLevel == 3) System.out.println("AlphaWeighted[I][2*L] ="+(int)AlphaWeighted[I][2*L]);
								if(VerboseLevel == 3) System.out.println("K = "+K);
								if((int)AlphaWeighted[I][2*L] == K){
									AlphaToUse=AlphaWeighted[I][2*L+1];
									Found=true;
									break; //exit L loop
								}//end if(AlphaWeighted[I][2*L] == K)
							}else{
								AlphaToUse = Alpha;
								Found = true;
								break;//exit L loop
							}//end if(Alpha == 0.5)
						}//end L loop
						if(Found == false){
							System.out.println("  Error -- alpha not found in Wiggle()");
							System.exit(1);
						}
						if(VerboseLevel == 3) System.out.println("AlphaToUse = "+AlphaToUse);
						Test = WC * E / ((1 - AlphaToUse) * V);
						if(VerboseLevel == 3) System.out.println("Length = "+Length+" Test = "+Test);
						if(Length > Test){
							if(VerboseLevel == 3) System.out.println("Time = "+Time);
							if(DynParm == true || Time <= DelT){//don't output multiple times
								java.io.File outfile;
								try {
									File file = new File("Wiggle.dng");
                     							if (Time < DelT){
                     								if(file.exists()){
                     									file.delete(); // remove residual file
                     								}
                     							}
                     							FileWriter fw = new FileWriter(file.getName(),true);
                     							BufferedWriter bw = new BufferedWriter(fw); 
                     							if (Time == DelT){
                     								if (DynParm == true){
                     									bw.write("For State Variable"+ISV+ " At time "+Time+" positivity/wiggle criterion violated between compartments "+I+" and "+K+" Maximum length is "+Test);
                     								}else{
                     									bw.write("Positivity/wiggle criterion violated at all timesteps between compartments"+I+"and "+K+ " Maximum length is "+ Test);
										}
									   bw.write("\r\n");//line feed
									}
									bw.close();
								} catch (IOException e){
									 e.printStackTrace();
								}								
							}//end if(DynParm == true || Time <= DelT)
						}// end if(Length > Test)
					}//end if (Q<0)
				}//end J loop
			}//end if
		}// end I loop
	}//end ISV loop
//KWL 3/15 change -- deleted VB code below ...
}//end Wiggle


 static void Newton(double x[],double XPREV[],float CumTime) {
 	 if(VerboseLevel > 1) System.out.println("I'm in Newton() and CumTime = "+CumTime);
 	 int I,J,K,BisectionMax,BigK;
	 boolean IEND,LambdaOK;
	 double TrialX,Lambda;
	 float f,fTrial;
	 xTrial = new double[NEQN+1];
	 NewtStep = new double[NEQN+1];
/*
	 THE NEWTON ALGORITHM IS X(K+1) = X(K) - STEP(K)
	STEP = SOLUTION TO LINEAR SYSTEM XJAC(X(K))*STEP = G(X(K))
	SC IS THE STOPPING CRITERION
	X() IS THE SET OF UNKNOWNS AT THE CURRENT TIME
*/
	try {
		File file = new File("NewtonInfo.dng");
		if(CumTime < DelT){
			if(file.exists()){
				file.delete();//remove residual file
			}
		}
		FileWriter fw = new FileWriter(file.getName(),true);
		BufferedWriter bw = new BufferedWriter(fw);	
		bw.write("\r\n");//line feed
		bw.write("Time = "+CumTime);
		bw.write("\r\n");//line feed
		bw.close();
	} catch (IOException e){
		 e.printStackTrace();
	}	

	IEND = false;
	K = 0;
	BigK = 0;
	do {
		K++;
		System.out.println("Newton iteration: "+K);
		if (K > MaxNewtonIteration){
			BigK +=1;
			if(BigK <= MaxSamplingIteration){
				if(SaveNewton == true){
					try {
						File file = new File("NewtonInfo.dng");
						FileWriter fw = new FileWriter(file.getName(),true);
						BufferedWriter bw = new BufferedWriter(fw);	
						bw.write("Random sampling iteration: "+BigK);
						bw.write("\r\n");//line feed
						bw.close();
					} catch (IOException e){
						 e.printStackTrace();
					}
				}
				//try best of randomly generated trial values
				//fTrialMin = Math.pow(1.0,30);
				double fTrialMin = 1.0e+30;
				for (I=1;I<NumSamples+1;I++){
					if(I == 1) {
						System.out.println("Randomizing ...");
						if(SaveNewton == true){
							try {
								File file = new File("NewtonInfo.dng");
								FileWriter fw = new FileWriter(file.getName(),true);
								BufferedWriter bw = new BufferedWriter(fw);	
								bw.write("Randomizing ...");
								bw.write("\r\n");//line feed
								bw.close();
							} catch (IOException e){
								e.printStackTrace();
							}
						}
					}
					for(J=1;J<NEQN+1;J++){
						xTrial[J] = LB+Math.random()*(UB-LB);
						if (VerboseLevel == 3) System.out.println("xTrial["+J+"] = "+xTrial[J]);
					}
					System.out.println("Calling EvalG");
					EvalG(xTrial,XPREV,G,999,A,BBC,b,CumTime);
					fTrial = 0;
					for(J=1;J<NEQN+1;J++){
						fTrial+=(float) Math.pow(G[J],2)/2;
					}
					if(fTrial<fTrialMin){
						fTrialMin=fTrial;
						for(J=1;J<NEQN+1;J++){
							x[J]=xTrial[J];
						}
					}
					f= (float)fTrialMin;
					K=1;//reset
				}// end I loop
				if(SaveNewton == true){
						try {
							File file = new File("NewtonInfo.dng");
							FileWriter fw = new FileWriter(file.getName(),true);
							BufferedWriter bw = new BufferedWriter(fw);	
							bw.write("Best f of random samples = "+fTrialMin);
							bw.write("\r\n");//line feed
							bw.close();
						} catch (IOException e){
							e.printStackTrace();
						}
					}
			}else{//BigK > MaxSamplingIteration
				if(SaveNewton == true){
					try {
						File file = new File("NewtonInfo.dng");
						FileWriter fw = new FileWriter(file.getName(),true);
						BufferedWriter bw = new BufferedWriter(fw);	
						bw.write("Algorithm terminated due to excessive Newton method iterations");
						bw.write("\r\n");//line feed
						bw.close();
					} catch (IOException e){
						e.printStackTrace();
					}
				}
				System.out.println("Algorithm terminated due to excessive Newton method iterations");
				System.exit(1);
			}
		}else{// K < MaxNewtonIteration
			EvalG(x,XPREV,G,999,A,BBC,b,CumTime);
			//check stopping criterion
			IEND = false;
			f=0;//initialize G(transpose)G quadratic function
			for(I=1;I<NEQN+1;I++){
				 f+=Math.pow(G[I],2)/2;
				 //System.out.println("G["+I+"] = "+G[I]);
			}
			
			//System.out.println("f = "+f+" and SC = "+SC);
			if(f<=SC){
				IEND = true;
			}
			//System.out.println("IEND = "+IEND);
			if (SaveNewton == true){
				try {
					File file = new File("NewtonInfo.dng");
					FileWriter fw = new FileWriter(file.getName(),true);
					BufferedWriter bw = new BufferedWriter(fw); 
					bw.write("Newton Iteration: "+K);
					bw.write("\r\n");//line feed
					bw.write("Current X values ...");
					bw.write("\r\n");//line feed
					for(I=1;I<NEQN+1;I++){
						String dummy = Double.toString(x[I]);
						bw.write(dummy);
						if(I<NEQN) bw.write(", ");
					}
					bw.write("\r\n");//line feed
					
					bw.write("Current f = "+f/NewtonSF);
					bw.write("\r\n");//line feed	
					bw.close();
				} catch (IOException e){
					 e.printStackTrace();
				}
			}
			    if(IEND == true) {
			    	    if(VerboseLevel == 3) System.out.println("IEND = true in Newton, exiting algorithm");
			    	    break;//exit do
			    }
			    
			    Jacobian(DeltFD,x,XPREV,G,XJAC,A,BBC,b,CumTime);
			    if (SaveNewton == true){
				try {
					File file = new File("NewtonInfo.dng");
					FileWriter fw = new FileWriter(file.getName(),true);
					BufferedWriter bw = new BufferedWriter(fw); 
					bw.write("Jacobian ...");
					bw.write("\r\n");//line feed

					for(I=1;I<NEQN+1;I++){
						for(J=1;J<NEQN+1;J++){
							String dummy = Double.toString(XJAC[I][J]);
							bw.write(dummy);
							if(J<NEQN) bw.write(", ");
						}
						bw.write("\r\n");//line feed
					}
					//bw.write("\r\n");//line feed
					bw.close();
				} catch (IOException e){
					 e.printStackTrace();
				}
;
			}

			//FIND THE NEWTON STEP -- RETURNED AS VECTOR NewtStep
			    for(I=1;I<NEQN+1;I++){
			    	    NewtStep[I] = G[I];
			    }
			  
			  if(VerboseLevel == 3){
			  	for(I=1;I<NEQN+1;I++){
			  		System.out.println("before LinearSolver NewtStep["+I+"] = "+NewtStep[I]);
			  		for (J=1;J<NEQN+1;J++){
			  			System.out.println("before LinearSolver XJAC["+I+"]["+J+"] = "+XJAC[I][J]);
			  		}
			  	}
			  }
			  // SOLVE LINEAR SYSTEM XJAC(X(K))*STEP = G(X(K))
			  LinearSolver(NEQN,XJAC,NewtStep);

			  // UPDATE X
			  Lambda = 1; // initialize at full Newton step
			  
			  boolean ExitLambdaLoop = false;
			  do{
			  	  //new! below!
			  	 LambdaOK = true;//initialize 
			  	  if(SaveNewton == true){
					  try {
						File file = new File("NewtonInfo.dng");
						FileWriter fw = new FileWriter(file.getName(),true);
						BufferedWriter bw = new BufferedWriter(fw); 
						bw.write("Lambda = " +Lambda);
						bw.write("\r\n");//line feed
						bw.close();
					} catch (IOException e){
						 e.printStackTrace();
					}
				  }
			  	for(I = 1;I<NEQN+1;I++){// check for NegProblem
			  		xTrial[I] = x[I]-Lambda*NewtStep[I];
			  		if(NegProblem == true && xTrial[I] <0){
			  			LambdaOK = false;
			  			break;//new!
			  		}
			  	}
			  	if (LambdaOK == true){
					//check for f improvement
					EvalG(xTrial,XPREV,G,999,A,BBC,b,CumTime);
					fTrial = 0;
					for(I=1;I<NEQN+1;I++){
						fTrial+=Math.pow(G[I],2)/2;
					}
					if(SaveNewton == true){
						try {
							File file = new File("NewtonInfo.dng");
							FileWriter fw = new FileWriter(file.getName(),true);
							BufferedWriter bw = new BufferedWriter(fw); 
							bw.write("Trial f = "+fTrial/NewtonSF);
							bw.write("\r\n");//line feed
							bw.close();
						} catch (IOException e){
							 e.printStackTrace();
						}
					}
					if(fTrial >= f)LambdaOK = false;
				}
				if(LambdaOK == false){
					Lambda = Lambda/2;
					if(Lambda < Epsilon){
						K = MaxNewtonIteration;
						ExitLambdaLoop = true;// exit Lambda loop
					}
				}else{
					for(I=1;I<NEQN+1;I++){
						x[I]=xTrial[I];
					}
					ExitLambdaLoop = true;
				}
			  }while(ExitLambdaLoop == false);//end ExitLambdaLoop loop
			  if(SaveNewton == true){
				  try {
					File file = new File("NewtonInfo.dng");
					FileWriter fw = new FileWriter(file.getName(),true);
					BufferedWriter bw = new BufferedWriter(fw); 
					bw.write("Updated X values ...");
					bw.write("\r\n");//line feed
					for (I=1;I<NEQN+1;I++){
						String dummy = Double.toString(x[I]);
						bw.write(dummy);
						if(I<NEQN) bw.write(", ");
					}
					bw.write("\r\n");//line feed
					bw.close();
				} catch (IOException e){
					 e.printStackTrace();
				}
			  }
 		}
	} while (IEND == false  );//end K loop
 }// end Newton()

 static void EvalG(double x[],double XPREV[],double G[],int ISELECTG,double A[][],double BBC[],double b[],float CumTime) {

 	if(VerboseLevel > 1) System.out.println("I'm in EvalG()"); 
 	 /*
	SUB TO EVALUATE G(X) = 0 FOR NEWTON'S METHOD
	xprev() is x at time 't'
	x() is x at time 't+1'
	*/

	int IEND,ISTART,I,J,ISELECT;
	float TempTime = 0;
	double ACFBtp1[];
	
	if(VerboseLevel == 3){
		for(I = 1;I<NEQN+1;I++){
			System.out.println("in EvalG, x[ "+I+"] = "+x[I]);
			
			System.out.println("in EvalG, XPREV[ "+I+"] = "+XPREV[I]);
		}
	}
	
	if(ISELECTG == 999) {//evaluate all G equations
		ISTART = 1;
		IEND = NEQN;
	}else{//evaluate only ISELECT equation
		ISTART = ISELECTG;
		IEND = ISELECTG;
	}
	
	if(VerboseLevel == 3)System.out.println("ISTART = "+ISTART+" and IEND = "+IEND);
	if(VerboseLevel == 3)System.out.println("NT = "+NT);
	if(Dynamic == true){
		if(DynType == 1){//Euler
		
		          //g(c) = [R(c)Vc}^t+1 - [R(c)V]^t - delT[Ac + f(c) + b]^t = 0	
			RVC = new double[NEQN+1];
			//compute RVC at time t
			
			if (CumTime == 0 || DynParm == true){
				ParmUpdate(CumTime);
			} 
			if(Rnonlinear == true) RofC(XPREV,CumTime);
			for(I=ISTART;I<IEND+1;I++){
			
				if(EquationSolver == false){
					RVC[I]=R[I]*V[I]*XPREV[I];// XPREV = time t
				}else{
				System.out.println("VB code not translated in EvalG()");
				//  RVC(I) = RVEqnSolver(I) * XPREV(I)
				}
			} 
			//compute ACFB at time t
			if(EquationSolver == false){
				if(CumTime == 0 || DynParm == true){
				
					//load transport terms into array A
					LoadATrans();//BBC loaded with boundary conditions only and on LHS
					//load linear source/sink terms
					LoadALinearSrcSnk();
				}
			}
			if(OutputTime == true){
				TempTime = CumTime;
				Aout(TempTime);// KWL 3/14 change
				RVout(TempTime);//KWL 3/14 change
			}

			//get loads and f at time t
			TempTime = CumTime;
			ForceFn(TempTime);
			FofC(XPREV,TempTime);
			if(VerboseLevel == 3){
				for(I = 1;I<NEQN+1;I++){
					System.out.println("in EvalG after FofC call, f["+I+"] = "+f[I]);
						}
			}
			if(Rnonlinear == true) RofC(XPREV,TempTime);
			
			ACFB = new double[NEQN+1];
			for(I=ISTART;I<IEND+1;I++){
				ACFB[I]=0;
				for(J=1;J<NEQN+1;J++){
					ACFB[I]+=A[I][J]*XPREV[J];
				}
				ACFB[I]+=f[I]+b[I];
			}
			
			RVCtp1 = new double[NEQN+1];
			//compute RVC at time t+1
			if (DynParm == true){
				ParmUpdate(CumTime+DelT);
			} 	
			TempTime = CumTime + DelT;
			if(Rnonlinear == true) RofC(x,TempTime);
			for(I=ISTART;I<IEND+1;I++){
				if(EquationSolver == false){
					RVCtp1[I]=R[I]*V[I]*x[I];// X = time t+1
				}else{
				System.out.println("VB code not translated in EvalG()");
				//  RVCtp1(I) = RVEqnSolver(I) * x(I)
				}
			}	
			//Evaluate G
			for(I=ISTART;I<IEND+1;I++){
				G[I]=RVCtp1[I]-RVC[I]-DelT*ACFB[I];
				if(VerboseLevel == 3)System.out.println("G["+I+"] = "+G[I]);
			}
		}// end Euler
		if(DynType == 2){
			System.out.println("MacCormack not appropriate for nonlinear system");
			System.exit(1);
		}
		if(DynType == 3){// BT
	
			//g(c) = [R(c)Vc]^t+1 - [R(c)V]^t - delT[Ac + f(c) + b]^t+1 = 0
            
			RVC = new double[NEQN+1];
			//compute RVC at time t
			if (CumTime == 0 || DynParm == true){
				ParmUpdate(CumTime);
			} 
			if(Rnonlinear == true) {
				TempTime = CumTime;
				RofC(XPREV,CumTime);
			}
			for(I=ISTART;I<IEND+1;I++){
				if(EquationSolver == false){
					RVC[I]=R[I]*V[I]*XPREV[I];// XPREV = time t
				}else{
				System.out.println("VB code not translated in EvalG()");
				//  RVC(I) = RVEqnSolver(I) * XPREV(I)
				}
			}
			RVCtp1 = new double[NEQN+1];
			//compute RVC at time t+1
			if (DynParm == true){
				ParmUpdate(CumTime+DelT);
				}
			if (Rnonlinear == true) {
				TempTime = CumTime+DelT;
				RofC(x,TempTime);
			}
			for(I=ISTART;I<IEND+1;I++){
				if(EquationSolver == false){
					RVCtp1[I]=R[I]*V[I]*x[I];// X = time t+1
				}else{
				System.out.println("VB code not translated in EvalG()");
				//  RVCtp1(I) = RVEqnSolver(I) * x(I)
				}
			}	
			//compute ACFB at time t+1
			if(EquationSolver == false){
				if(CumTime == 0 || DynParm == true){
					//load transport terms into array A
					LoadATrans();//BBC loaded with boundary conditions only and on LHS
					//load linear source/sink terms
					LoadALinearSrcSnk();
				}
			}
		
			if(OutputTime == true){
				TempTime = CumTime;
				Aout(TempTime);// KWL 3/14 change
				RVout(TempTime);//KWL 3/14 change
			}
			//get loads and f at time t+1
			TempTime = CumTime+DelT;
			ForceFn(TempTime);
			FofC(x,TempTime);
			if(VerboseLevel == 3){
				for(I = 1;I<NEQN+1;I++){
					System.out.println("in EvalG after FofC call, f["+I+"] = "+f[I]);
						}
			}
			ACFBtp1 = new double[NEQN+1];
			for(I=ISTART;I<IEND+1;I++){
				ACFBtp1[I]=0;
				for(J=1;J<NEQN+1;J++){
					ACFBtp1[I]+=A[I][J]*x[J];
				}
				ACFBtp1[I]+=f[I]+b[I];
			}
			//Evaluate G
			for(I=ISTART;I<IEND+1;I++){
				G[I]=RVCtp1[I]-RVC[I]-DelT*ACFBtp1[I];
				if(VerboseLevel == 3)System.out.println("G["+I+"] = "+G[I]);
			}
		}// end BT method
		if(DynType == 4){//CT
			//g(c) = [R(c)Vc]^t+1 - [R(c)V]^t - (delT/2){[Ac + f(c) + b]^t+1 + [Ac + f(c) + b]^t}= 0	
			RVC = new double[NEQN+1];//R*V*C at time t
			ACFB = new double[NEQN+1];//A*c+f+b at time t
			
			//compute RVC and ACFB at time t
			if (CumTime == 0 || DynParm == true){
				ParmUpdate(CumTime);
			} 	
			if(Rnonlinear == true) {
				RofC(XPREV,CumTime);
				}
			for(I=ISTART;I<IEND+1;I++){
				if(EquationSolver == false){
					RVC[I]=R[I]*V[I]*XPREV[I];// XPREV = time t
				}else{
					System.out.println("VB code not translated in EvalG()");
				//  RVC(I) = RVEqnSolver(I) * XPREV(I)
				}
			}
			//get loads and f at 't'
			ForceFn(CumTime);
			FofC(XPREV,CumTime);
			ACFB = new double[NEQN+1];
			for(I=ISTART;I<IEND+1;I++){
				ACFB[I]=0;
				for(J=1;J<NEQN+1;J++){
					ACFB[I]+=A[I][J]*XPREV[J];
				}
				ACFB[I]+=f[I]+b[I];
			}
			
			RVCtp1 = new double[NEQN+1];//R*V*C at time t+1
			ACFBtp1 = new double[NEQN+1];//A*c+f+b at time t+1
			//compute RVC and ACFB at time t+1
			if (CumTime == 0 || DynParm == true){
				ParmUpdate(CumTime+DelT);
			} 	
			if(Rnonlinear == true) RofC(x,CumTime+DelT);
			for(I=ISTART;I<IEND+1;I++){
				if(EquationSolver == false){
					RVCtp1[I]=R[I]*V[I]*x[I];// x = time t+1
				}else{
					System.out.println("VB code not translated in EvalG()");
				//   RVCtp1(I) = RVEqnSolver(I) * x(I)
				}
			}
			if(EquationSolver == false){
				if(CumTime == 0 || DynParm == true){
					//load transport terms into array A
					LoadATrans();//BBC loaded with boundary conditions only and on LHS
					//load linear source/sink terms
					LoadALinearSrcSnk();
				}
			}
		
			if(OutputTime == true){
				TempTime = CumTime;
				Aout(TempTime);// KWL 3/14 change
				RVout(TempTime);//KWL 3/14 change
			}
						
			//get loads and f at 't+1'
			TempTime = CumTime+DelT;
			ForceFn(TempTime);
			FofC(x,TempTime);
			

			for(I=ISTART;I<IEND+1;I++){
				ACFBtp1[I]=0;
				for(J=1;J<NEQN+1;J++){
					ACFBtp1[I]+=A[I][J]*x[J];
				}
				ACFBtp1[I]+=f[I]+b[I];
			}
			
			//Evaluate G
			for(I=ISTART;I<IEND+1;I++){
				G[I]=RVCtp1[I]-RVC[I]-(DelT/2)*(ACFBtp1[I]+ACFB[I]);
				if(VerboseLevel == 3)System.out.println("G["+I+"] = "+G[I]);
			}
		}//end CT method
		
	}else{//steady state
		FofC(x,0);
		//Evaluate G
		for(I=ISTART;I<IEND+1;I++){
			G[I] = 0;
			for(J=1;J<NEQN+1;J++){
				G[I]+=A[I][J]*x[J];
			}
			G[I] += f[I]-b[I];
			if(VerboseLevel == 3)System.out.println("G["+I+"] = "+G[I]);
		}
	}//end steady state
 }//end EvalG()
 
 static void Jacobian(double DeltFD,double x[],double XPREV[],double G[],double XJAC[][],double A[][],double BBC[],double b[],float CumTime){
 	if(VerboseLevel == 3) System.out.println("I'm in Jacobian()");
 	int I,J,K,Col;
 	double Diff1,Diff2,XT[],GT[];
 	XT = new double[NEQN+1];
 	GT = new double[NEQN+1];
	double dfdc[][];
	dfdc = new double[NEQN+1][NEQN+1];	
	double dRdc[];
	dRdc = new double[NEQN+1];

 
 	if(AnalDeriv == false){
 
 		for(I=1;I<NEQN+1;I++){
 			XT[I]=x[I];//initialize
 		}
 		
 		for(I=1;I<NEQN+1;I++){
		System.out.println("Estimating row "+I+" of Jacobian matrix ...");	
 			for(J=1;J<NEQN+1;J++){
 				// constructing dG[I]/dx[J]
 				XT[J]=x[J]+DeltFD*x[J];// forward difference
 				//System.out.println("XT["+J+"] = "+XT[J]);
 				EvalG(XT,XPREV,GT,I,A,BBC,b,CumTime);
 				Diff1 = GT[I];
 				XT[J]=x[J]-DeltFD*x[J];// backward difference
 				EvalG(XT,XPREV,GT,I,A,BBC,b,CumTime);
 				Diff2 = GT[I];
 				XJAC[I][J]=(Diff1-Diff2)/(2*DeltFD*x[J]); // central difference
 				//System.out.println("in Jacobian, XJAC["+I+"]["+J+"] = "+XJAC[I][J]);
 			}//end J loop
 		}// end I loop
  	} else{//AnalDeriv == true
 		Getdfdc(dfdc,x,CumTime);
 		if(Dynamic == true){
 			GetdRdC(dRdc, x,CumTime);
 			if(DynType == 1){//Euler
 				for(I=1;I<NEQN+1;I++){
 					Col = 0;
 					for(J=1;J<MaxNSV+1;J++){//walk state variables
 						for(K=1;K<NT+1;K++){//walk interior compartments
 							if (Compartment[K][3] == 0 && SVComp[J][K] == 1){// interior compartment
 								Col++;
 								if(I==Col){
 									XJAC[I][Col]=V[Col]*(R[Col]+x[I]*dRdc[Col]);
 								}else{
 									XJAC[I][Col] = 0;
 								}
 							}
 						}
 					}
 				}
 			}
 			if(DynType == 2){//MacCormack
 				System.out.println("Error: MacCormack not implemented for nonlinear case");
 				System.exit(1);
 			}
 			if(DynType == 3 || DynType == 4){//BT and CT
 				for(I=1;I<NEQN+1;I++){
 					Col = 0;
 					for(J=1;J<MaxNSV+1;J++){//walk state variables
 						for(K=1;K<NT+1;K++){//walk interior compartments
 							if (Compartment[K][3] == 0 && SVComp[J][K] == 1){// interior compartment
 								Col++;
 								if(I==Col){
 									if(DynType == 3) XJAC[I][Col]=V[Col]*(R[Col]+x[I]*dRdc[Col])-DelT*(A[I][Col]+dfdc[I][Col]);
 									if(DynType == 4) XJAC[I][Col]=V[Col]*(R[Col]+x[I]*dRdc[Col])-(DelT/2)*(A[I][Col]+dfdc[I][Col]);
 								}else{
 									if(DynType == 3) XJAC[I][Col]=-DelT*(A[I][Col]+dfdc[I][Col]);
 									if(DynType == 4) XJAC[I][Col]=-(DelT/2)*(A[I][Col]+dfdc[I][Col]);
 								}
 							}
 						}
 					}
 				}
	
 			}
 		}else{//steady state
 			for(I=1;I<NEQN+1;I++){
 				for(J=1;J<NEQN+1;J++){
 					XJAC[I][J]=A[I][J]+dfdc[I][J];
 				}
 			}
 		}
 	}//end if(AnalDeriv == false)
 	//System.out.println("Leaving Jacobian");
  }//end Jacobian()
  
  static void TransportFluxToBoundaryCompartments(float TempTime,int ISV,int AdjComp,double Flux){
		if(VerboseLevel > 1) System.out.println("I'm in TransportFluxToBoundaryCompartments()");
		int I,J;
 
		try {
                    File file = new File("TransportFluxToBoundaryCompartments.csv");
                    
                   if (TempTime < DelT){ 
                         if(file.exists()){
                        file.delete(); // remove residual file
                        }
                    }
                    FileWriter fw = new FileWriter(file.getName(),true);
                    BufferedWriter bw = new BufferedWriter(fw); 
                    if (TempTime < DelT){
                        bw.write("Time, SV, BoundaryCompartment, Flux");
			bw.write("\r\n");//line feed
                    }
                    if(TempTime >= DelT) {
                    	bw.write(TempTime+","+ISV+","+AdjComp+","+Flux);    
                    	bw.write("\r\n");//line feed
                    }
                 bw.close();
		} catch (IOException e){
			 e.printStackTrace();
		}	
	}// end TransportFluxToBoundaryCompartments method
	
	static void ArealFluxToBoundaryCompartments(float TempTime,int ISV,int AdjComp,double Flux){
		if(VerboseLevel > 1) System.out.println("I'm in ArealFluxToBoundaryCompartments()");
		int I,J;

		try {
                    File file = new File("ArealFluxToBoundaryCompartments.csv");
                    
                   if (TempTime < DelT){ 
                         if(file.exists()){
                        file.delete(); // remove residual file
                        }
                    }
                    FileWriter fw = new FileWriter(file.getName(),true);
                    BufferedWriter bw = new BufferedWriter(fw); 
                    if (TempTime < DelT){
                        bw.write("Time, SV, BoundaryCompartment, Flux");
 			bw.write("\r\n");//line feed
                    }
                    if(TempTime >= DelT) {
                    	bw.write(TempTime+","+ISV+","+AdjComp+","+Flux);    
                    	bw.write("\r\n");//line feed
                    }
                 bw.close();
		} catch (IOException e){
			System.out.println("Error: Exception in ArealFluxToBoundaryCompartment() method");
			 e.printStackTrace();
			 
		}	
	}// end ArealFluxToBoundaryCompartments method

 static void FofC(double x[],double TempTime) {
	if(VerboseLevel > 1) System.out.println("I'm in FofC()");
	int I;
	for (I=1;I<NEQN+1;I++){
		f[I] = 0;
	}
	//application-specific code
	/* following code for nonlinear kinetics example 4.4 in GEM book
	double Kgmax,Ks;
	Kgmax = 5;
	Ks = 150;
	f[1] = Kgmax*x[2]*V[2]*x[1]/(Ks+x[2]);
	f[2]=-f[1];
	//end application-specific code
	if(VerboseLevel == 3){
		for(I=1;I<NEQN+1;I++){
			System.out.println("in FofC, x["+I+"] ="+x[I]);
			System.out.println("in FofC, f["+I+"] ="+f[I]);
		}
	}
	//*/
	//end application-specific code
	
	if(Shellfofc == true) {
		if(Dynamic == true){
			WriteTimeAndStatus("fofc", TempTime, x);
		}else{
			WriteStatus("fofc", x);
		}
		ShellFunction("fofc.exe");
		
		try {
		       String strFile = "fofc.csv";
		       BufferedReader br = new BufferedReader( new FileReader(strFile));//create BufferedReader to read csv file
		       String dummystring = "";
		       StringTokenizer st = null;
		       if(VerboseLevel > 1)System.out.println("Reading fofc.csv");
		       for (I=1; I < NEQN+1;I++) {
				 dummystring = br.readLine();
				 f[I]=Double.parseDouble(dummystring);
				 if(VerboseLevel > 1) System.out.println("f["+I+"] = "+f[I]);
			 } // end for i
			 br.close(); // close input stream
		}// end try
		catch(Exception e)
		 {
			System.out.println("Exception while reading fofc.csv file: " + e);
		 }// end catch	
		 
		 if(OutputTime == true){
		 	 try {
				FileWriter fw = new FileWriter("FofC.dng");
				BufferedWriter out = new BufferedWriter(fw);
				out.write("At time "+TempTime);
				out.write("\r\n");// line feed
				out.write("For following C values ...");
				out.write("\r\n");// line feed
				for (I = 1; I < NEQN+1;I++){
					String dummy = Double.toString(x[I]);
					out.write(dummy);
				out.write("\r\n");// line feed
				}
				out.write("result in following F values ...");
				out.write("\r\n");// line feed
				for (I = 1; I < NEQN+1;I++){
					String dummy = Double.toString(f[I]);
					out.write(dummy);
				out.write("\r\n");// line feed
				}
				out.close();
			 } catch (IOException e){
			 	 e.printStackTrace();
			 }
		 }
	} // end if(Shellfofc == true)

 } // end fofC method	
 
 	static void Getdfdc(double dfdc[][],double x[],double TempTime){
		if (VerboseLevel > 1) System.out.println("I'm in GetdfdC");
		int I,J;
		for(I=1;I<NEQN+1;I++){
			for(J=1;J<NEQN+1;J++){
				dfdc[I][J]=0;//initialize
			}
		}
		if(Shellfofc == true){
			if(Dynamic == true){
				WriteTimeAndStatus("dfdc", TempTime, x);
			}else{
				WriteStatus("dfdc", x);
			}
			if(EquationSolver == true){
				System.out.println("VB code not translated in Getdfdc method");
			/*
				Call ShellFunction("dfdc.exe")
					Open "dfdc.csv" For Input As #1
						For I = 1 To NEQN
							For J = 1 To NEQN
								Input #1, dfdc(I, J)
							Next J
						Next I
					Close #1
			*/
			}else{
				ShellFunction("dfdc.exe");
				try {
				       String strFile = "dfdc.csv";
				       BufferedReader br = new BufferedReader( new FileReader(strFile));//create BufferedReader to read csv file
				       String dummystring = "";
				       StringTokenizer st = null;
				       if(VerboseLevel > 1)System.out.println("Reading dfdc.csv");
				       for (I=1; I < NEQN+1;I++) {
				       	       for(J=1;J<NEQN+1;J++){
						 dummystring = br.readLine();
						 dfdc[I][J]=Double.parseDouble(dummystring);
						 if(VerboseLevel > 1) System.out.println("dfdc["+I+"]["+J+"]= "+dfdc[I][J]);
					       }//end for J
					 } // end for i
					 br.close(); // close input stream
				}// end try
				catch(Exception e)
				 {
					System.out.println("Exception while reading dfdc.csv file: " + e);
				 }// end catch	
        		}//end if(EquationSolver == true)

		}//end if(Shellfofc == true)
		
		//application-specific code
		/* following code for nonlinear kinetics example 4.4 in GEM book
	
		double Kgmax,Ks;
		Kgmax = 5;
		Ks = 150;
		dfdc[1][1]= Kgmax*x[2]*V[2]/(Ks+x[2]);
		System.out.println("dfdc11 = "+dfdc[1][1]);
		dfdc[1][2] = -x[1]*Kgmax*x[2]*V[2]/(Math.pow((Ks+x[2]),2))+(Kgmax*V[2]*x[1])/(Ks+x[2]);
		dfdc[2][1]=-dfdc[1][1];
		dfdc[2][2]=-dfdc[1][2];
		*/
		//end application-specific code
	
	}//end GetDFdC method
	
	static void RofC(double x[],double TempTime){
		if(VerboseLevel > 1) System.out.println("I'm in RofC()");
		int I;
		for(I = 1;I<NEQN+1;I++){
			R[I] = 0;//initialize
		}
		if(ShellRofc == true){
			WriteTimeAndStatus("Rofc", TempTime, x);
			ShellFunction("Rofc.exe");
			try {
			       String strFile = "Rofc.csv";
			       BufferedReader br = new BufferedReader( new FileReader(strFile));//create BufferedReader to read csv file
			       String dummystring = "";
			       StringTokenizer st = null;
			       if(VerboseLevel > 1)System.out.println("Reading fofc.csv");
			       for (I=1; I < NEQN+1;I++) {
					 dummystring = br.readLine();
					 R[I]=Double.parseDouble(dummystring);
					 if(VerboseLevel > 1) System.out.println("R["+I+"] = "+R[I]);
				 } // end for i
				 br.close(); // close input stream
			}// end try
			catch(Exception e)
			{
				System.out.println("Exception while reading fofc.csv file: " + e);
			}// end catch	

			 if(OutputTime == true){
				 try {
					FileWriter fw = new FileWriter("RofC.dng");
					BufferedWriter out = new BufferedWriter(fw);
					out.write("At time "+TempTime);
					out.write("\r\n");// line feed
					out.write("For following C values ...");
					out.write("\r\n");// line feed
					for (I = 1; I < NEQN+1;I++){
						String dummy = Double.toString(x[I]);
						out.write(dummy);
					out.write("\r\n");// line feed
					}
					out.write("result in following R values ...");
					out.write("\r\n");// line feed
					for (I = 1; I < NEQN+1;I++){
						String dummy = Double.toString(R[I]);
						out.write(dummy);
					out.write("\r\n");// line feed
					}
					out.close();
				 } catch (IOException e){
					 e.printStackTrace();
				 }
			 }
		}//end if(ShellRofc == true)
		
		//application-specific code
		/*
		//following code for nonlinear sorption example 4.5 in GEM book
		double BulkDensity = 2500000;
		double WaterContent = 0.5;
		double KdNonlinear = 0.000000135;
		double KdExponent = 0.5;
		int I;
		for(I=1;I<NEQN+1;I++){
			R[I]=1+(BulkDensity*KdNonlinear/WaterContent)*Math.pow(x[I],(KdExponent-1));
		}
		//end application-specific code
		*/
	}//end RofC method
	
	static void GetdRdC(double dRdC[],double x[],double TempTime){
		if(VerboseLevel > 1) System.out.println("I'm in GetdRdC()");
		int I;
		for(I=1;I<NEQN+1;I++){
			dRdC[I]=0;//initialize
		}
		if(ShellRofc == true){
			WriteTimeAndStatus("dRdc", TempTime, x);

			if(EquationSolver == true){
				System.out.println("VB code not translated in Getdfdc method");
			/*
				Call ShellFunction("dRdc.exe")
					Open "dfdc.csv" For Input As #1
						For I = 1 To NEQN
							Input #1, dRdc(I)
						Next I
					Close #1
			*/
			}else{
				ShellFunction("dRdc.exe");
				try {
				       String strFile = "dRdc.csv";
				       BufferedReader br = new BufferedReader( new FileReader(strFile));//create BufferedReader to read csv file
				       String dummystring = "";
				       StringTokenizer st = null;
				       if(VerboseLevel > 1)System.out.println("Reading dfdc.csv");
				       for (I=1; I < NEQN+1;I++) {
						 dummystring = br.readLine();
						 dRdC[I]=Double.parseDouble(dummystring);
						 if(VerboseLevel > 1) System.out.println("dRdc["+I+"]= "+dRdC[I]);
					 } // end for i
					 br.close(); // close input stream
				}// end try
				catch(Exception e)
				 {
					System.out.println("Exception while reading dRdc.csv file: " + e);
				 }// end catch	
        		}//end if(EquationSolver == true)

		}//end if(Shellfofc == true)

		//application-specific code
		/*
		double BulkDensity = 2500000;
		double WaterContent = 0.5;
		double KdNonlinear = 0.000000135;
		double KdExponent = 1.5;
		int I;

		for(I=1;I<NEQN+1;I++){
			dRdC[I]=(BulkDensity*KdNonlinear/WaterContent)*(KdExponent-1)*Math.pow(x[I],(KdExponent-2));
		}
		*/
		//end application-specific code
	}// end GetdRdC method
	

	
	static void ShellFunction(String FileToFire){
		if(VerboseLevel >1) System.out.println("I'm in ShellFunction()and FileToFire = "+FileToFire);
      	       Runtime run = Runtime.getRuntime();
      	       Process pr = null;
      	       try {
      	       	       pr = run.exec(FileToFire);
      	       } catch (IOException e) {
      	       	       e.printStackTrace();
      	       }
      	       //now wait for process to run
      	       try{
      	       	       pr.waitFor();
      	       } catch(InterruptedException e){
      	       	       System.out.println("Error: exception in ShellFunction() method");
      	       	       e.printStackTrace();
      	       }
	} // end ShellFunction method
	
	static void WriteTimeAndStatus(String filename, Double Time, double c[]){
		if(VerboseLevel >1) System.out.println("I'm in WriteTimeAndStatus and filename = "+filename);
		try {
			FileWriter fw = new FileWriter("Shell"+filename+".csv");
			BufferedWriter out = new BufferedWriter(fw);
			String dummy = Double.toString(Time);
			out.write(dummy);
			out.write("\r\n");// line feed
			for (int I = 1; I < NEQN+1;I++){
				dummy = Double.toString(c[I]);
				out.write(dummy);
				out.write("\r\n");// line feed
			}
			out.close();
		} catch (IOException e){
			 System.out.println("Error: exception in WriteTimeAndStatus() method");
			 e.printStackTrace();
		}
	} //end WriteTimeAndStatus

	static void WriteStatus(String filename,double x[]){
		if(VerboseLevel >1) System.out.println("I'm in WriteStatusand filename = "+filename);
		try {
			FileWriter fw = new FileWriter("Shell"+filename+".csv");
			BufferedWriter out = new BufferedWriter(fw);
			for (int I = 1; I < NEQN+1;I++){
				String dummy = Double.toString(x[I]);
				out.write(dummy);
				out.write("\r\n");// line feed
			}
			out.close();
		} catch (IOException e){
			 System.out.println("Error: exception in WriteStatus() method");
			 e.printStackTrace();
		}
	} //end WriteStatus
	
	static void WriteTime(String filename, double Time){
		if(VerboseLevel >1) System.out.println("I'm in WriteTimes and filename = "+filename);
		try {
			FileWriter fw = new FileWriter("Shell"+filename+".csv");
			BufferedWriter out = new BufferedWriter(fw);
			String dummy = Double.toString(Time);
			out.write(dummy);
			out.write("\r\n");// line feed
			out.close();
		} catch (IOException e){
			 System.out.println("Error: exception in WriteTime() method");
			 e.printStackTrace();
		}
	} //end WriteTimeAndStatus
	
}// end Class
