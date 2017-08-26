package fun.grn.grneat.evaluators;

import fun.grn.grneat.evolver.GRNGenome;
import fun.grn.grneat.grn.GRNModel;

public class LowPassFilterEvaluator extends GRNGenomeEvaluator {
	double coef=3.0;

	public LowPassFilterEvaluator() {
		name = "MichalSignalProcessExp3";
		
		numGRNInputs=1;
		numGRNOutputs=1;
		
	}

	@Override
	public double evaluate(GRNGenome aGenome) {
		GRNModel grn=buildGRNFromGenome(aGenome);
		double fitness=evolveGRNWithThreshold(grn, false, 125.0, 1000)+
				evolveGRNWithThreshold(grn, false, 250.0, 1000)+
				evolveGRNWithThreshold(grn, false, 75, 1000)+
				evolveGRNWithThreshold(grn, false, 31.25, 1000)+
				evolveGRNDoubleFreq(grn, false, 125, 31.25, 1000)+
				evolveGRNDoubleFreq(grn, false, 250, 31.25, 1000)+
				evolveGRNZero(grn, false, 1000);

		//System.err.println("fitness="+fitness+"  =>  "+fun.grn.grneat.grn.toString());
		aGenome.setNewFitness(fitness);

		GRNGenomeEvaluator.numEvaluations++;
		return fitness;
	}

	double evolveGRNZero(GRNModel grn, boolean printTrace, int nStepMax) {
		double fitness=0.0;

		grn.reset();
		grn.evolve(25);
		for (int nStep=0; nStep<nStepMax; nStep++) {
			double dt;
			grn.proteins.get(0).concentration=0.0;
			grn.evolve(1);
			dt=0;
			fitness+=Math.abs(dt-Math.min(1.0,grn.proteins.get(1).concentration*coef))*(1+dt);
			if (printTrace) {
				System.out.println(nStep+"\t"+Math.min(1.0,grn.proteins.get(1).concentration*coef)+"\t"+dt+"\t"+0.0);
			}
			if (grn.proteins.get(1).concentration!=grn.proteins.get(1).concentration) {
				System.err.println("nan");
			}
		}

		//System.err.println("fitness="+fitness+"  =>  "+fun.grn.grneat.grn.toString());
		return -fitness;

	}

	double evolveGRNWithThreshold(GRNModel grn, boolean printTrace, double halfFreq, int nStepMax) {
		double fitness=0.0;

		grn.reset();
		grn.evolve(25);
		double ots[]={0,0,0,0,0,0,0,0,0,0};
		int nbEvents=0;
		int lastEvent=0;
		for (int nStep=0; nStep<nStepMax; nStep++) {
			double dt;
			double in=Math.sin(nStep*Math.PI/halfFreq-Math.PI/2)/2+0.5;
			grn.proteins.get(0).concentration=in;
			grn.evolve(1);
			if (halfFreq<=100) {
				dt=0;
			} else {
				dt=Math.sin(nStep*Math.PI/halfFreq-Math.PI/2)/2+0.5;;
			}
			for (int nst=0; nst<9; nst++) {
				ots[nst]=ots[nst+1];
			}
			ots[9]=Math.min(1.0,grn.proteins.get(1).concentration*coef);
			if (nStep-lastEvent>10 && ((ots[0]<0.5 && ots[9]>=0.5) || (ots[0]>=0.5 && ots[9]<0.5))) {
				nbEvents++;
				//System.out.println("Event "+lastEvent+"  "+nStep);
				lastEvent=nStep;
			}
			fitness+=Math.abs(dt-Math.min(1.0,ots[9]))*(1+dt);
			if (printTrace) {
				System.out.println(nStep+"\t"+Math.min(1.0,ots[9])+"\t"+dt+"\t"+in);
			}
			if (ots[9]!=ots[9]) {
				System.err.println("nan");
			}
		}
		double S;
		double nbEventDesired=halfFreq<=100?0:(int)(2*(double)nStepMax/halfFreq);
		if (nbEvents==0 || nbEvents>nbEventDesired*2) {
			S=0;
		} else {
			S=1.0-(double)(Math.abs(nbEvents-nbEventDesired))/nbEventDesired;
		}
		if (printTrace) {
//			System.out.println(nbEvents+"  "+ S+ "  " +nbEventDesired);
		}
		fitness*=1.0/(1.0+S);

		//System.err.println("fitness="+fitness+"  =>  "+fun.grn.grneat.grn.toString());
		return -fitness;
	}

	double evolveGRNDoubleFreq(GRNModel grn, boolean printTrace, double halfFreq1, double halfFreq2, int nStepMax) {
		double fitness=0.0;

		grn.reset();
		grn.evolve(25);
		double ots[]={0,0,0,0,0,0,0,0,0,0};
		int nbEvents=0;
		int lastEvent=0;
		for (int nStep=0; nStep<nStepMax; nStep++) {
			double dt;
			double in=0.25*Math.sin(nStep*Math.PI/halfFreq1-Math.PI/2)+0.25*Math.sin(nStep*Math.PI/halfFreq2-Math.PI/2)+0.5;
			grn.proteins.get(0).concentration=in;
			grn.evolve(1);
			for (int nst=0; nst<9; nst++) {
				ots[nst]=ots[nst+1];
			}
			dt=0.5*Math.sin(nStep*Math.PI/halfFreq1-Math.PI/2)+0.5;;
			ots[9]=Math.min(1.0,grn.proteins.get(1).concentration*coef);
			if (nStep-lastEvent>10 && ((ots[0]<0.5 && ots[9]>=0.5) || (ots[0]>=0.5 && ots[9]<0.5))) {
				nbEvents++;
				//System.out.println("Event "+lastEvent+"  "+nStep);
				lastEvent=nStep;
			}
			fitness+=Math.abs(dt-Math.min(1.0,ots[9]))*(1+dt);
			if (printTrace) {
				System.out.println(nStep+"\t"+Math.min(1.0,ots[9])+"\t"+dt+"\t"+in);
			}
			if (ots[9]!=ots[9]) {
				System.err.println("nan");
			}
		}
		double S;
		double nbEventDesired= (int)(2*(double)nStepMax/halfFreq1);
		if (nbEvents==0 || nbEvents>nbEventDesired*2) {
			S=0;
		} else {
			S=1.0-(double)(Math.abs(nbEvents-nbEventDesired))/nbEventDesired;
		}
		if (printTrace) {
//			System.out.println(nbEvents+"  "+ S+ "  " +nbEventDesired);
		}
		fitness*=1.0/(1.0+S);

		//System.err.println("fitness="+fitness+"  =>  "+fun.grn.grneat.grn.toString());
		return -fitness;
	}

	double evolveGRNTwoFrequencies(GRNModel grn, boolean printTrace, double halfFreq1, double halfFreq2, int nStepMax, int startF2, int endF2) {
		double fitness=0.0;

		grn.reset();
		grn.evolve(25);
		double ots[]={0,0,0,0,0,0,0,0,0,0};
		int nbEvents=0;
		int lastEvent=0;
		for (int nStep=0; nStep<nStepMax; nStep++) {
			double dt,in;
			if (nStep<startF2 || nStep>endF2) {
				in=Math.sin(nStep*Math.PI/halfFreq1-Math.PI/2)/2+0.5;
				dt=in;
			} else {
				in=Math.sin(nStep*Math.PI/halfFreq2-Math.PI/2)/2+0.5;
				dt=0;
			}
			grn.proteins.get(0).concentration=in;
			grn.evolve(1);
			for (int nst=0; nst<9; nst++) {
				ots[nst]=ots[nst+1];
			}
			ots[9]=Math.min(1.0,grn.proteins.get(1).concentration*coef);
			if (nStep-lastEvent>10 && ((ots[0]<0.5 && ots[9]>=0.5) || (ots[0]>=0.5 && ots[9]<0.5))) {
				nbEvents++;
				//System.out.println("Event "+lastEvent+"  "+nStep);
				lastEvent=nStep;
			}
			fitness+=Math.abs(dt-Math.min(1.0,ots[9]))*(1+dt);
			if (printTrace) {
				System.out.println(nStep+"\t"+Math.min(1.0,ots[9])+"\t"+dt+"\t"+in);
			}
			if (ots[9]!=ots[9]) {
				System.err.println("nan");
			}
		}
		double S;
		double nbEventDesired= (int)(2*(double)nStepMax/halfFreq1);
		if (nbEvents==0 || nbEvents>nbEventDesired*2) {
			S=0;
		} else {
			S=1.0-(double)(Math.abs(nbEvents-nbEventDesired))/nbEventDesired;
		}
		if (printTrace) {
//			System.out.println(nbEvents+"  "+ S+ "  " +nbEventDesired);
		}
		fitness*=1.0/(1.0+S);

		//System.err.println("fitness="+fitness+"  =>  "+fun.grn.grneat.grn.toString());
		return -fitness;
	}

	public static void main(String args[]) throws Exception {
		LowPassFilterEvaluator eval=new LowPassFilterEvaluator();
		double freqs1[]=	{25, 	50, 	66,		75};
		double freqs2[]=	{125, 	200, 	500, 	750};
		int nSteps1[]=	{1000, 	1000, 	1000, 	1000};
		int nSteps2[]=	{1000, 	1000, 	2000, 	5000};
		GRNModel greatGRN[]=new GRNModel[26];
		GRNModel gaGRN[]=new GRNModel[26];
		GRNModel esGRN[]=new GRNModel[26];
		String greatFiles[]={
				"grn_299_-329.09981349860095.fun.grn.grneat.grn",	"grn_299_-406.21324950142565.fun.grn.grneat.grn",	"grn_299_-465.3027738438956.fun.grn.grneat.grn",
				"grn_299_-329.2325839380702.fun.grn.grneat.grn",	"grn_299_-409.1461016744854.fun.grn.grneat.grn",	"grn_299_-482.22024647740636.fun.grn.grneat.grn",
				"grn_299_-346.2301485332818.fun.grn.grneat.grn",	"grn_299_-409.36748107116995.fun.grn.grneat.grn",	"grn_299_-509.73837272483564.fun.grn.grneat.grn",
				"grn_299_-347.7666559154336.fun.grn.grneat.grn",	"grn_299_-411.2289911450258.fun.grn.grneat.grn",	"grn_299_-511.6519835240472.fun.grn.grneat.grn",
				"grn_299_-348.77350451429635.fun.grn.grneat.grn",	"grn_299_-425.2215394406491.fun.grn.grneat.grn",	"grn_299_-513.2001161028146.fun.grn.grneat.grn",
				"grn_299_-369.5373808236058.fun.grn.grneat.grn",	"grn_299_-428.96734035951727.fun.grn.grneat.grn",	"grn_299_-531.584740500517.fun.grn.grneat.grn",
				"grn_299_-369.96784305136396.fun.grn.grneat.grn",	"grn_299_-434.43962622292213.fun.grn.grneat.grn",	"grn_299_-574.6193158328872.fun.grn.grneat.grn",
				"grn_299_-388.77514852472444.fun.grn.grneat.grn",	"grn_299_-436.38405426543727.fun.grn.grneat.grn",	"grn_299_-585.3843633836484.fun.grn.grneat.grn",
				"grn_299_-397.7053490688321.fun.grn.grneat.grn",	"grn_299_-453.29152264146956.fun.grn.grneat.grn",	"grn_299_-628.0302373539503.fun.grn.grneat.grn",
				"grn_299_-405.68737821830524.fun.grn.grneat.grn",	"grn_299_-456.3915636292151.fun.grn.grneat.grn",	"grn_299_-984.4006489695837.fun.grn.grneat.grn",
		};
		String gaFiles[]={
				"grn_300_-1092.3221329658606.fun.grn.grneat.grn",	"grn_300_-533.2285004737496.fun.grn.grneat.grn",	"grn_300_-707.681543613565.fun.grn.grneat.grn",
				"grn_300_-439.34021964103385.fun.grn.grneat.grn",	"grn_300_-543.0755120436231.fun.grn.grneat.grn",	"grn_300_-721.0833806784334.fun.grn.grneat.grn",
				"grn_300_-460.49726960695455.fun.grn.grneat.grn",	"grn_300_-547.0213301001882.fun.grn.grneat.grn",	"grn_300_-755.9050783342359.fun.grn.grneat.grn",
				"grn_300_-495.26687479449987.fun.grn.grneat.grn",	"grn_300_-549.8234310668693.fun.grn.grneat.grn",	"grn_300_-812.9594210821596.fun.grn.grneat.grn",
				"grn_300_-497.4030481033237.fun.grn.grneat.grn",	"grn_300_-563.8182109096965.fun.grn.grneat.grn",	"grn_300_-837.171211773604.fun.grn.grneat.grn",
				"grn_300_-500.517446421063.fun.grn.grneat.grn",	"grn_300_-564.9190505161973.fun.grn.grneat.grn",	"grn_300_-862.7678862624423.fun.grn.grneat.grn",
				"grn_300_-502.4635133364742.fun.grn.grneat.grn",	"grn_300_-588.2281512496638.fun.grn.grneat.grn",	"grn_300_-886.1584632306841.fun.grn.grneat.grn",
				"grn_300_-506.0805448594242.fun.grn.grneat.grn",	"grn_300_-630.5825397450902.fun.grn.grneat.grn",	"grn_300_-942.707677192221.fun.grn.grneat.grn",
				"grn_300_-512.0730671589673.fun.grn.grneat.grn",	"grn_300_-664.2426636730385.fun.grn.grneat.grn",	"grn_300_-979.7518431256908.fun.grn.grneat.grn",
				"grn_300_-515.4858299560933.fun.grn.grneat.grn",	"grn_300_-687.3888845088879.fun.grn.grneat.grn",	"grn_300_-989.7703156656927.fun.grn.grneat.grn",
		};
		String esFiles[]={
				"grn_299_-1048.3973279031354.fun.grn.grneat.grn",	"grn_299_-376.18310431612105.fun.grn.grneat.grn",	"grn_299_-433.70418871061565.fun.grn.grneat.grn",
				"grn_299_-1078.6455970297197.fun.grn.grneat.grn",	"grn_299_-376.4382630147122.fun.grn.grneat.grn",	"grn_299_-441.73550612575474.fun.grn.grneat.grn",
				"grn_299_-1105.173333125703.fun.grn.grneat.grn",	"grn_299_-390.0920998088891.fun.grn.grneat.grn",	"grn_299_-442.11965421818604.fun.grn.grneat.grn",
				"grn_299_-332.6502728189408.fun.grn.grneat.grn",	"grn_299_-402.6369767692752.fun.grn.grneat.grn",	"grn_299_-454.5296590729756.fun.grn.grneat.grn",
				"grn_299_-337.9259564761536.fun.grn.grneat.grn",	"grn_299_-415.23823410693836.fun.grn.grneat.grn",	"grn_299_-486.5677766751558.fun.grn.grneat.grn",
				"grn_299_-345.5411007646542.fun.grn.grneat.grn",	"grn_299_-419.3784431378082.fun.grn.grneat.grn",	"grn_299_-502.1929000488882.fun.grn.grneat.grn",
				"grn_299_-350.4317540971891.fun.grn.grneat.grn",	"grn_299_-419.69150481850005.fun.grn.grneat.grn",	"grn_299_-521.6252954327341.fun.grn.grneat.grn",
				"grn_299_-351.0796597899462.fun.grn.grneat.grn",	"grn_299_-422.7547119276268.fun.grn.grneat.grn",	"grn_299_-709.0543870179614.fun.grn.grneat.grn",
				"grn_299_-361.52229849842996.fun.grn.grneat.grn",	"grn_299_-426.4821064514153.fun.grn.grneat.grn",	"grn_299_-727.3027958269574.fun.grn.grneat.grn",
				"grn_299_-368.2405556367121.fun.grn.grneat.grn",	"grn_299_-426.9370165895936.fun.grn.grneat.grn",	"grn_299_-965.8748184524359.fun.grn.grneat.grn",
		};
		for (int i=0; i<greatGRN.length; i++) {
			greatGRN[i]=GRNModel.loadFromFile("/Users/cussat/Recherche/Projets/grnNEAT/GREAT_GIT/launcher016-hyperion006/Generalization/LPF/GREAT/"+greatFiles[i]);
			gaGRN[i]=GRNModel.loadFromFile("/Users/cussat/Recherche/Projets/grnNEAT/GREAT_GIT/launcher016-hyperion006/Generalization/LPF/GA/"+gaFiles[i]);
			esGRN[i]=GRNModel.loadFromFile("/Users/cussat/Recherche/Projets/grnNEAT/GREAT_GIT/launcher016-hyperion006/Generalization/LPF/ES/"+esFiles[i]);
		}
		
		
		System.out.println("GREAT");
		for (int j=0; j<greatGRN.length; j++) {
			
			System.out.print("\t"+greatFiles[j]);
		}
		System.out.println();
		for (int i=0; i<freqs1.length; i++) {
			System.out.print(freqs1[i]);
			for (int j=0; j<greatGRN.length; j++) {
				double greatCurFit=-eval.evolveGRNWithThreshold(greatGRN[j], false, freqs1[i], nSteps1[i]);
				System.out.print("\t"+greatCurFit);
			}
			System.out.println();
		}
		for (int i=0; i<freqs2.length; i++) {
			System.out.print(freqs2[i]);
			for (int j=0; j<greatGRN.length; j++) {
				double greatCurFit=-eval.evolveGRNWithThreshold(greatGRN[j], false, freqs2[i], nSteps2[i]);
				System.out.print("\t"+greatCurFit);
			}
			System.out.println();
		}
		for (int i=0; i<freqs2.length; i++) {
			for (int j=0; j<freqs1.length; j++) {
				System.out.print(freqs2[i]+" ("+freqs1[j]+") ");
				for (int k=0; k<greatGRN.length; k++) {
					double greatCurFit=-eval.evolveGRNDoubleFreq(greatGRN[k], false, freqs2[i], freqs1[i], nSteps2[i]);
					System.out.print("\t"+greatCurFit);
				}
				System.out.println();
			}
		}
/*		for (int i=0; i<freqs2.length; i++) {
			for (int j=0; j<freqs1.length; j++) {
				double greatFit=0;
				double gaFit=0;
				System.out.print(freqs2[i]+"_"+freqs1[j]+"_"+freqs2[i]);
				for (int k=0; k<greatGRN.length; k++) {
					double greatCurFit=-eval.evolveGRNTwoFrequencies(greatGRN[k], false, freqs2[i], freqs1[i], nSteps2[i], 1000, 2000);
					double gaCurFit=-eval.evolveGRNTwoFrequencies(gaGRN[k], false, freqs2[i], freqs1[i], nSteps2[i], 1000, 2000);
					greatFit+=greatCurFit;
					gaFit+=gaCurFit;
					System.out.print("\t"+greatCurFit+"\t"+gaCurFit);
					greatAvg[k]+=greatCurFit;
					gaAvg[k]+=gaCurFit;
				}
				System.out.println("\t"+greatFit/greatGRN.length+"\t"+gaFit/greatGRN.length);
			}
		}*/
		
		
		System.out.println("GA");
		for (int j=0; j<gaGRN.length; j++) {
			
			System.out.print("\t"+gaFiles[j]);
		}
		System.out.println();
		for (int i=0; i<freqs1.length; i++) {
			System.out.print(freqs1[i]);
			for (int j=0; j<gaGRN.length; j++) {
				double gaCurFit=-eval.evolveGRNWithThreshold(gaGRN[j], false, freqs1[i], nSteps1[i]);
				System.out.print("\t"+gaCurFit);
			}
			System.out.println();
		}
		for (int i=0; i<freqs2.length; i++) {
			System.out.print(freqs2[i]);
			for (int j=0; j<gaGRN.length; j++) {
				double gaCurFit=-eval.evolveGRNWithThreshold(gaGRN[j], false, freqs2[i], nSteps2[i]);
				System.out.print("\t"+gaCurFit);
			}
			System.out.println();
		}
		for (int i=0; i<freqs2.length; i++) {
			for (int j=0; j<freqs1.length; j++) {
				System.out.print(freqs2[i]+" ("+freqs1[j]+") ");
				for (int k=0; k<gaGRN.length; k++) {
					double gaCurFit=-eval.evolveGRNDoubleFreq(gaGRN[k], false, freqs2[i], freqs1[i], nSteps2[i]);
					System.out.print("\t"+gaCurFit);
				}
				System.out.println();
			}
		}
		
		
		System.out.println("ES");
		for (int j=0; j<esGRN.length; j++) {
			
			System.out.print("\t"+esFiles[j]);
		}
		System.out.println();
		for (int i=0; i<freqs1.length; i++) {
			System.out.print(freqs1[i]);
			for (int j=0; j<esGRN.length; j++) {
				double esCurFit=-eval.evolveGRNWithThreshold(esGRN[j], false, freqs1[i], nSteps1[i]);
				System.out.print("\t"+esCurFit);
			}
			System.out.println();
		}
		for (int i=0; i<freqs2.length; i++) {
			System.out.print(freqs2[i]);
			for (int j=0; j<esGRN.length; j++) {
				double esCurFit=-eval.evolveGRNWithThreshold(esGRN[j], false, freqs2[i], nSteps2[i]);
				System.out.print("\t"+esCurFit);
			}
			System.out.println();
		}
		for (int i=0; i<freqs2.length; i++) {
			for (int j=0; j<freqs1.length; j++) {
				System.out.print(freqs2[i]+" ("+freqs1[j]+") ");
				for (int k=0; k<esGRN.length; k++) {
					double esCurFit=-eval.evolveGRNDoubleFreq(esGRN[k], false, freqs2[i], freqs1[i], nSteps2[i]);
					System.out.print("\t"+esCurFit);
				}
				System.out.println();
			}
		}
				
		System.exit(0);		
		
/*		String grnFile="default.fun.grn.grneat.grn";
		double halfFreq1=125;
		double halfFreq2=31.25;
		int steps=1000;
		int mode=0;
		
		for (int k=0; k<args.length; k++) {
			if (args[k].compareTo("mode") == 0) {
				k++;
				mode=Integer.parseInt(args[k]);
			} else if (args[k].compareTo("fun.grn.grneat.grn")==0) {
				k++;
				grnFile=args[k];
			} else if (args[k].compareTo("halfFreq1")==0) {
				k++;
				halfFreq1=Double.parseDouble(args[k]);
			} else if (args[k].compareTo("halfFreq2")==0) {
				k++;
				halfFreq2=Double.parseDouble(args[k]);
			} else if (args[k].compareTo("step")==0) {
				k++;
				steps=Integer.parseInt(args[k]);
			} else {
				System.out.println("Unrecognized option: "+args[k]);
			}
		}
		GRNModel fun.grn.grneat.grn = GRNModel.loadFromFile(grnFile);
		MichalSignalProcessExp3 eval=new MichalSignalProcessExp3();
		double fit=0;
		switch (mode) {
		case 0:
			fit=eval.evolveGRNWithThreshold(fun.grn.grneat.grn, true, halfFreq1, steps);
			break;
		case 1:
			fit=eval.evolveGRNDoubleFreq(fun.grn.grneat.grn, true, halfFreq1, halfFreq2, steps);
			break;
		case 2:
			fit=eval.evolveGRNTwoFrequencies(fun.grn.grneat.grn, true, halfFreq1, halfFreq2, steps, 1000, 2000);
			break;
		case 3:
			fit=eval.evolveGRNZero(fun.grn.grneat.grn, true, steps);
		default:
			System.out.println("Unrecognized mode: "+mode);
			break;
		}
		
//		double fit=eval.evolveGRNWithThreshold(fun.grn.grneat.grn, true, 125.0, 1000);
//		double fit=eval.evolveGRNWithThreshold(fun.grn.grneat.grn, true, 250.0, 1000);
//		double fit=eval.evolveGRNWithThreshold(fun.grn.grneat.grn, true, 75, 1000);
//		double fit=eval.evolveGRNWithThreshold(fun.grn.grneat.grn, true, 31.25, 1000);
//		double fit=eval.evolveGRNDoubleFreq(fun.grn.grneat.grn, true, 125, 31.25, 1000);
//		double fit=eval.evolveGRNDoubleFreq(fun.grn.grneat.grn, true, 250, 31.25, 1000);
//		double fit=eval.evolveGRNZero(fun.grn.grneat.grn, true, 1000);
//		System.out.println("Fitness="+fit);*/
	}

}
