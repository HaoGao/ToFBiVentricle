#include "tetVolumeCalculation.h" 


double tetVolumeCalculation(tetgenio& out)
{
     //output the information for debugging
	int nElem = out.numberoftetrahedra;
	int nPoint = out.numberofpoints;
	
	//std::cout<<"nPoint: "<<nPoint<<"\n";
	//std::cout<<"nElem:  "<<nElem<<"\n";
	
	
	REAL *pointlist = out.pointlist;
	int *tetrahedronlist = out.tetrahedronlist;
	
	//print out the pointlist
	/*std::cout<<"output the nodes\n";
	int pstartIndex = 0;
	for (int pindex = 0; pindex< nPoint; pindex++)
	{
	    pstartIndex = pindex*3;
		std::cout<<pointlist[pstartIndex]<<"\t";
		std::cout<<pointlist[pstartIndex+1]<<"\t";
		std::cout<<pointlist[pstartIndex+2]<<"\n";
	}*/
	
	//print out the elem list
	//std::cout<<"\n\n output the elements\n";
	double volTotal = 0.0;
	int elemStartIndex = 0;
	for (int elemIndex = 0; elemIndex < nElem; elemIndex ++)
	{
		elemStartIndex = elemIndex*4; //only support first order tet
		/*std::cout << tetrahedronlist[elemStartIndex]<<"\t";
		std::cout << tetrahedronlist[elemStartIndex+1]<<"\t";
		std::cout << tetrahedronlist[elemStartIndex+2]<<"\t";
		std::cout << tetrahedronlist[elemStartIndex+3]<<"\n";
		*/
		unsigned int n1, n2, n3, n4;
		n1 = tetrahedronlist[elemStartIndex];
		n2 = tetrahedronlist[elemStartIndex+1];
		n3 = tetrahedronlist[elemStartIndex+2];
		n4 = tetrahedronlist[elemStartIndex+3];
		
		//extract the coordinates
		double p1[3], p2[3], p3[3], p4[3];
		
		//p1 node
		p1[0]= pointlist[n1*3];
		p1[1]= pointlist[n1*3+1];
		p1[2]= pointlist[n1*3+2];
		
		//p2 node
		p2[0]= pointlist[n2*3];
		p2[1]= pointlist[n2*3+1];
		p2[2]= pointlist[n2*3+2];
		
		//p3 node
		p3[0]= pointlist[n3*3];
		p3[1]= pointlist[n3*3+1];
		p3[2]= pointlist[n3*3+2];
		
		//p4 node
		p4[0]= pointlist[n4*3];
		p4[1]= pointlist[n4*3+1];
		p4[2]= pointlist[n4*3+2];
		
		//calculate the volume
		REAL a[3], u[3], v[3];
		a[0] = p1[0]-p4[0];
		a[1] = p1[1]-p4[1];
		a[2] = p1[2]-p4[2];
		
		u[0] = p2[0]-p4[0];
		u[1] = p2[1]-p4[1];
		u[2] = p2[2]-p4[2];
		
		v[0] = p3[0]-p4[0];
		v[1] = p3[1]-p4[1];
		v[2] = p3[2]-p4[2];
		
		double uxv[3];
		uxv[0] =  u[1]*v[2] - u[2]*v[1];
        uxv[1] = -u[0]*v[2] + u[2]*v[0];
        uxv[2] =  u[0]*v[1] - u[1]*v[0];
		
	    double vol = 0.0;
	    vol = a[0]*uxv[0] + a[1]*uxv[1] + a[2]*uxv[2];
        volTotal = volTotal + fabs(vol)/6.0;
		
	}
	
	
	
	return volTotal;
}


double tetVolumeCalculation(std::vector< std::vector<double> >& endo_points, int NoOfEndoNode)
{
	tetgenio in, out;
	 tetgenbehavior b;
	 
	 in.numberofpoints = NoOfEndoNode;
     in.pointlist = new double[in.numberofpoints * 3];
   
       for (unsigned int pIndex = 0; pIndex<NoOfEndoNode; pIndex++)
      {				
		in.pointlist[pIndex*3+0] = endo_points[pIndex][0];
	    in.pointlist[pIndex*3+1] = endo_points[pIndex][1];
	    in.pointlist[pIndex*3+2] = endo_points[pIndex][2];	
	  }	
	
     char switches[] = {"Q"};
  	 if (!b.parse_commandline(switches)) {
        terminatetetgen(NULL, 10);
       }
	
	
	//in.save_nodes("endo");
	tetrahedralize(&b, &in, &out);
	double volTotal = 0.0;
	volTotal = tetVolumeCalculation(out);
	
	return volTotal;
}
