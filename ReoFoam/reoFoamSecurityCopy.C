/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    pimpleFoam

Description
    Transient solver for incompressible, turbulent flow of Newtonian fluids,
    with optional mesh motion and mesh topology changes.

    Turbulence modelling is generic, i.e. laminar, RAS or LES may be selected.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "singlePhaseTransportModel.H"
#include "kinematicMomentumTransportModel.H"
#include "pimpleControl.H"
//#include "CorrectPhi.H"
#include "fvOptions.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"
    #include "createDyMControls.H"
    #include "createFields.H"
    #include "createUfIfPresent.H"

  //turbulence->validate();

    if (!LTS)
    {
        #include "CourantNo.H"
        #include "setInitialDeltaT.H"
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        #include "readDyMControls.H"
      // correctPhi=true;
        if (LTS)
        {
            #include "setRDeltaT.H"
        }
        else
        {
            #include "CourantNo.H"
            #include "setDeltaT.H"
        }


	mesh.update();

	//	MRF.update();

	//phi = mesh.Sf() & Uf();

	//	fvc::makeRelative(phi, U);

	if (checkMeshCourantNo)
	      {
#include "meshCourantNo.H"
	      }

	//	   gamma.correctBoundaryConditions();
 
	  // Vector
	  //	  vector ComponenteY(0, 1, 0);
	//	 MRF.correctBoundaryVelocity(U);

	  // dimensionedSymmTensor  ComponenteY =  fvc::grad(fvc::curl(U)));
	  
			      //	  fvOptions.constrain(Uprel);
			      //	  Uprel.solve();
	  //		  B =fvc::curl(U);

	//     gamma.correctBoundaryConditions();	

	     forAll(U, celli){

	     	      U[celli].z() = 0.0;

	     	    }
	     //   fvc::correctUf(Uf, U, phi);
	     phi =  mesh.Sf() & Uf();
	   fvc::makeRelative(phi, U);

	//	correctUf(Uf,U, phi);
	//	phi = fvc::absolute(phi,U)
	    fvScalarMatrix gammaEqn
	      (
	       //fvm::ddt(gamma) +fvc::div(phi ,gamma) // -beta*fvm::laplacian(gamma)+beta*fvc::laplacian(gamma)
	       fvm::ddt(gamma)  + fvc::div(phi, gamma) //-fvm::Sp(fvc::div(fvc::meshPhi(U)),gamma)  
	       // fvm::ddt(gamma) + gamma*fvc::div(U) + fvc::Sp(1.0,U)()&fvc::div(gamma)()
	       );
	    //	    gammaEqn.setReference(0, 1.0);	    
	    // gammaEqn.relax();
	    //   fvOptions.constrain(gammaEqn);

	    //  gammaTransport = fvc::div(phi,gamma);
	    //    fvc::makeAbsolute(phi, U);
	    //
	tmp<volTensorField> tgradU(fvc::grad(U));
	//	const volTensorField& gradU = tgradU();

	const volTensorField& gradU = tgradU();

	


	//	vector lul(1.0, 1.0, 0.0);
	//volSymmTensorField Dshear2(diagT(gradU));
	//		volTensorField Dext = D*inv(D).ref();
	//	gradG = fvc::grad(gamma);				      //	Dext.value().xx=0.0;
	//term1 = Theta*fvc::curl(fvc::curl(U)());
	//term2 = fvc::div(NodiagT(gradU));
	//    phi = mesh.Sf() & Uf();

	
	  //dRau2.primitiveFieldRef() = (Tensor)diag(dRau2.primitiveFieldRef()).ref();
	  // Dshear.replace(0,0.0);
	  // Dshear.replace(4,0.0);
	  // Dshear.replace(8,0.0);
	  /////
	  //
	//	  Dext = Dext-Dshear;
	    //    phi = fvc::interpolate(U)() ;  
	//	MRF.correctBoundaryVelocity(U);
	// pimple.correctNonOrthogonal();

		//	    fvc::correctUf(Uf, U, phi);
	     fvVectorMatrix UEqn
	    (
	     //	     fvm::laplacian(U)==-Ma*fvc::grad(gamma) +  fvc::div(fvc::curl(U))
	     // 


	  
	     //Theta*fvc::grad(fvc::div(U)())
	    
	     //	      2*fvc::div(diagT(gradU))
		      // -fvc::grad(fvc::div(U)()))


	       //  -fvc::laplacian(U)-fvc::curl(fvc::curl(U)()))

	       // 	        Ma*fvc::grad(gamma)+
	       //   Theta*(fvm::laplacian(U)+fvc::curl(fvc::curl(U)()))
	       //   +Tr/4.0*(fvc::div(diagT(gradU))-fvc::grad(fvc::div(U)))
	       // +2*fvc::div(NodiagT(gradU))


	     //Buena
	     fvm::laplacian(Theta+Tr/4.0, U)// - beta*fvm::ddt(U)+beta*fvc::ddt(U)
	     //==-Ma*fvc::grad(gamma)  
	     	   //  - 2.0*fvc::div(NodiagT(gradU))/(Theta+Tr/4.0)
	     //	     	     		      -Theta*fvc::curl(fvc::curl(U)())
	     //	   (Tr/4.0-1.0)*(2.0*fvc::div(NodiagT(gradU)))


	       // 	     (Tr/4.0)*fvm::laplacian(U)- beta*fvm::ddt(U)+beta*fvc::ddt(U)
	       //  ==-Ma*fvc::grad(gamma)  
 	       //     - 2.0*fvc::div(NodiagT(gradU))-Theta*fvc::grad(fvc::div(U)())
	       // +Tr/4.0*(2.0*fvc::div(NodiagT(gradU)))
		     
	     );
	     // gammaTransport = -Ma*fvc::grad(gamma);
	// UEqn.setReference(0,lul);
	  //	  fvVectorMatrix& UEqn = tUEqn.ref();
	    //	   	  	  UEqn.relax();
	    //	    UEqn.setReference(0, 1.0);

	    //    fvOptions.constrain(UEqn);

	     //	      UEqn.relax();

	     //UEqn.relax();

	     // UEqn.relax();
	     
	     forAll(U, celli){

	       U[celli].z() = 0.0;

	     }
	     
	         gammaEqn.relax();
	     gammaEqn.solve();	
	     UEqn.solve();

	     //	     fvc::makeRelative(phi, U);
	     // //	     fvc::makeRelative(phi, U);
	     // fvScalarMatrix gammaEqn2
	     //   (
	     // 	fvm::ddt(gamma)  + fvm::div(phi,gamma) //-fvm::Sp(fvc::div(fvc::meshPhi(U)),gamma)  
	     // 	);
	     // // fvc::makeAbsolute(phi, U);

	     // gammaEqn2.relax();
	     // gammaEqn2.solve();	
	     // UEqn.solve();

	     // //	     fvc::makeRelative(phi, U);
	     // fvScalarMatrix gammaEqn3
	     //   (
	     // 	fvm::ddt(gamma)  + fvm::div(phi,gamma) //-fvm::Sp(fvc::div(fvc::meshPhi(U)),gamma)  
	     // 	);
	     // //fvc::makeAbsolute(phi, U);
 	     // UEqn.solve();
	     // gammaEqn3.relax();
	     // gammaEqn3.solve();	
	     // UEqn.solve();

	     // fvc::makeRelative(phi, U);
	     // fvScalarMatrix gammaEqn4
	     //   (
	     // 	fvm::ddt(gamma)  + fvc::div(phi,gamma) //-fvm::Sp(fvc::div(fvc::meshPhi(U)),gamma)  
	     // 	);
	     // fvc::makeAbsolute(phi, U);
 	     // UEqn.solve();
	     // gammaEqn4.relax();
	     // gammaEqn4.solve();
	     


	     // 	     fvc::makeRelative(phi, U);
	     // fvScalarMatrix gammaEqn5
	     //   (
	     // 	fvm::ddt(gamma)  + fvc::div(phi,gamma) //-fvm::Sp(fvc::div(fvc::meshPhi(U)),gamma)  
	     // 	);
	     // fvc::makeAbsolute(phi, U);
 	     // UEqn.solve();
	     // gammaEqn5.relax();
	     // gammaEqn5.solve();	


	     // fvc::makeRelative(phi, U);
	     // fvScalarMatrix gammaEqn6
	     //   (
	     // 	fvm::ddt(gamma)  + fvc::div(phi,gamma) //-fvm::Sp(fvc::div(fvc::meshPhi(U)),gamma)  
	     // 	);
	     // fvc::makeAbsolute(phi, U);
 	     // UEqn.solve();
	     // gammaEqn6.relax();
	     // gammaEqn6.solve();


	     // fvc::makeRelative(phi, U);
	     // fvScalarMatrix gammaEqn7
	     //   (
	     // 	fvm::ddt(gamma)  + fvc::div(phi,gamma) //-fvm::Sp(fvc::div(fvc::meshPhi(U)),gamma)  
	     // 	);
	     // fvc::makeAbsolute(phi, U);
 	     // UEqn.solve();
	     // gammaEqn7.relax();
	     // gammaEqn7.solve();

	     // fvc::makeRelative(phi, U);
	     // fvScalarMatrix gammaEqn8
	     //   (
	     // 	fvm::ddt(gamma)  + fvc::div(phi,gamma) //-fvm::Sp(fvc::div(fvc::meshPhi(U)),gamma)  
	     // 	);
	     // fvc::makeAbsolute(phi, U);
 	     // UEqn.solve();
	     // gammaEqn8.relax();
	     // gammaEqn8.solve();	

	     //             gammaEqn.relax();

 
		    //		    gamma.correctBoundaryConditions();

	     // U.relax();

		    //   gammaEqn.relax()

	     //     gamma.relax();
	    // fvOptions.correct(gamma);
		   

	     //	fvOptions.constrain(UEqn);

	      //   fvOptions.correct(UEqn);

		  

	      //            UEqn.solve();

	     //U.correctBoundaryConditions();

	    //MRF.correctBoundaryVelocity(U);    // 
	    //	    fvc::corr ectUf(Uf, U, phi);


	    //	    phi =  mesh.Sf() & fvc::interpolate(U);
	    //phi =  fvc::interpolate(U);
	    // fvc::makeRelative(phi, U);
	    


	    //	    fvc::correctUf(Uf, U, phi);

	    // Make the fluxes relative to the mesh motion
	    //    fvc::makeRelative(phi, U);

	    
	    
	    // fvOptions.correct(gamma);
	    // surfaceScalarField sf = fvc::interpolate(Ma)& mesh.Sf();
	    // scalar Ma2 = 4.0;
	    // fvc::makeAbsolute(phi, U);

	    
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;


      Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
      label wallPatch = mesh.boundaryMesh().findPatchID("cylinder");
      
        runTime.write();

	if (runTime.write()){

	  //  angleV = beta*gamma;
	//   volVectorField firstEig = U;fvc::grad(gamma)() ;
	//   try {
	//     volTensorField T = fvc::grad(U)();eigenVectors(gradU);
	//     volTensorField Tu = fvc::grad(U)();eigenVectors(gradU);
	//     	firstEig.component(vector::X) = T.component(tensor::XX);
	//     auto& eig_ = U;
	//     scalar txx(1.0);
	//     forAll(T, celli) {
	//       //   eig_[celli]=eigenValues(Psi_[celli]);
	//       try {
		
	//   	Tu[celli][8] = txx;
	//   	eig_[celli] = eigenValues(Tu[celli]);
	//   	T[celli] = eigenVectors(Tu[celli], eig_[celli]);

	//       } 	catch(...)
	//   	{
	//   	}

	//     }
	//     forAll(firstEig, cellI)
	//       {

	//   		firstEig[cellI][0] = T[cellI][3];// 3  4  5
	//   	firstEig[cellI][1] = T[cellI][4];
	// 	firstEig[cellI][2] = T[cellI][5];
	//   	// if(mag(eig_[1])>mag(eig_[0]))
	//   	//   {
	//   	//     firstEig[cellI][0] = T[cellI][1];3  4  5
	//   	//     firstEig[cellI][1] = T[cellI][4];
	//   	//     firstEig[cellI][2] = T[cellI][7];
	//   	//   }
	//   	// else{
	//   	//     firstEig[cellI][0] = T[cellI][0];3  4  5
	//   	//     firstEig[cellI][1] = T[cellI][3];
	//   	//     firstEig[cellI][2] = T[cellI][6];
	//   	// }
	//   	// [0];
	//   	 Info<<"autovalores"<< eig_[cellI][0] <<"                "<< eig_[cellI][1]<<"                "<< eig_[cellI][2]<<"\n"; 
	//   	       Info<<firstEig[cellI][0] <<"                "<< firstEig[cellI][1]<<"                "<< firstEig[cellI][2]<<"\n"; 
		 
	//       }
	    

	    
	//     forAll(angleV,cellI)
	//       {
	//   	angleV[cellI] = Foam::sqrt(firstEig[cellI][0]*firstEig[cellI][0]+firstEig[cellI][1]*firstEig[cellI][1]);
	//   			Info<<angleV[cellI];

	//   	if(abs(angleV[cellI])>1e-8)
	//   	  {
	// 	    angleV[cellI] = firstEig[cellI][1]/angleV[cellI] *180/3.14159;//Pasamos a deg
	//   	  }
	//   	else{
	//   	  angleV[cellI]=0.0;
	//   	}
	//       }
	//   }
	//   catch(...)
	//     {
	//       forAll(angleV,cellI)
	//   	{
	//   	  		  angleV[cellI] = 0.0;//T[cellI].XY/angleV[cellI] *180/3.14159;//Pasamos a deg
	//   	}
	//     }

	  
	// angleV.write();
	  IOField<vector> cfOut
	    (
	     IOobject
	     (
	      "cf",
		mesh.time().timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::NO_WRITE
		),
	       mesh.Cf().boundaryField()[wallPatch]
	       // Sf ()normal   Cf center  magSf
	   );

	cfOut.write();


	IOField<tensor> gradOut
	  (
	       IOobject
	       (
		"gradU",
		mesh.time().timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::NO_WRITE
		),
	       gradU.boundaryField()[wallPatch]
	       // Sf ()normal   Cf center  magSf
	   );

	gradOut.write();




	
	  IOField<vector> SfOut
	  (
	       IOobject
	       (
		"sf",
		mesh.time().timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::NO_WRITE
		),
	       mesh.Sf().boundaryField()[wallPatch]
	       // Sf ()normal   Cf center  magSf
	   );

	  SfOut.write();

	  IOField<scalar> magSf
	  (
	       IOobject
	       (
		"magsf",
		mesh.time().timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::NO_WRITE
		),
	       mesh.magSf().boundaryField()[wallPatch]
	       // Sf ()normal   Cf center  magSf
	   );


	 magSf.write();

	 IOField<scalar> magGamma
	  (
	       IOobject
	       (
		"magGamma",
		mesh.time().timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::NO_WRITE
		),
	       gamma.boundaryField()[wallPatch]
	       // Sf ()normal   Cf center  magSf
	   );


	 magGamma.write();


	 
	
    }
	

    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
