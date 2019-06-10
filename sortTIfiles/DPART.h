/* MigSelect 2013-2019 Arun Sethuraman, Vitor Sousa, Jody Hey */
/* This code was developed based on IMa2  2009-2010  Jody Hey, Rasmus Nielsen and Sang Chul Choi*/
/* calculation of the mean partition (proposed by Huelsenbeck et al. 2007) ------------------------------------------------------------*/

   /*Cost matrix is generated using eq.16 in Konovalov et al. 2005 and solved by the hungarian algorithm.
   Original C++ source code of the hungarian algorithm is provided by 
           http://www.topcoder.com/tc?module=Static&d1=tutorials&d2=hungarianAlgorithm 
   Here it is converted to C */

#define Nv 50            /* max number of vertices in one part */
#define INF 100000000   /* just infinity */
#define whichmax(a, b) ((a) > (b) ? (a) : (b))
#define whichmin(a, b) ((a) < (b) ? (a) : (b))

//int cost[Nv][Nv];          /* cost matrix */
//int n, max_match;        /* n workers and n jobs */
//int lx[Nv], ly[Nv];        /* labels of X and Y parts */
//int xy[Nv];               /* xy[x] - vertex that is matched with x, */
//int yx[Nv];               /* yx[y] - vertex that is matched with y */
//int Sh[Nv], Th[Nv];        /* sets Sh and Th in algorithm */
//int slack[Nv];            /* as in the algorithm description */
//int slackx[Nv];           /* slackx[y] such a vertex, that */
//                         /* l(slackx[y]) + l(y) - w(slackx[y],y) = slack[y] */
//int prev[Nv];             /* array for memorizing alternating paths */

/*----------------------------------------------------------------------------------------------------*/
void add_to_tree(int x, int prevx); 
///* x - current vertex,prevx - vertex from X before x in the alternating path, */
///* so we add edges (prevx, xy[x]), (xy[x], x) */
//{
//	int y;
//
//    Sh[x] = 1;                      /* add x to Sh */
//    prev[x] = prevx;                /* we need this when augmenting */
//    for (y = 0; y < n; y++)         /* update slacks, because we add new vertex to Sh */
//        if (lx[x] + ly[y] - cost[x][y] < slack[y])
//        {
//            slack[y] = lx[x] + ly[y] - cost[x][y];
//            slackx[y] = x;
//        }
//}

/*----------------------------------------------------------------------------------------------------*/
void augment(void)  ;                       /* main function of the algorithm */
//{
//
//    int x, y, root;                        /* just counters and root vertex */
//    int q[Nv], wr = 0, rd = 0;              /* q - queue for bfs, wr,rd - write and read */
//	int cx,cy,ty;
//	int i, j;
//	int delta = INF;
//
//    if (max_match == n) return;            /* check wether matching is already perfect */
//                                           /* pos in queue */
//    memset(Sh, 0, sizeof(Sh));             /* init set Sh */
//    memset(Th, 0, sizeof(Th));             /* init set Th */
//    memset(prev, -1, sizeof(prev));        /* init set prev - for the alternating tree */
//
//
//    for (x = 0; x < n; x++)                /* finding root of the tree */
//        if (xy[x] == -1)
//        {
//            q[wr++] = root = x;
//            prev[x] = -2;
//            Sh[x] = 1;
//            break;
//        }
//
//    for (y = 0; y < n; y++)                /* initializing slack array */
//    {
//        slack[y] = lx[root] + ly[y] - cost[root][y];
//        slackx[y] = root;
//    }
//
//    /* second part of augment() function */
//    while (1)                                                           /* main cycle */
//    {
//        while (rd < wr)                                                 /* building tree with bfs cycle */
//        {
//            x = q[rd++];                                                /* current vertex from X part */
//            for (y = 0; y < n; y++)                                     /* iterate through all edges in equality graph */
//                if (cost[x][y] == lx[x] + ly[y] &&  !Th[y])
//                {
//                    if (yx[y] == -1) break;                             /* an exposed vertex in Y found, so */
//                                                                        /* augmenting path exists! */
//                    Th[y] = 1;                                          /* else just add y to Th, */
//                    q[wr++] = yx[y];                                    /* add vertex yx[y], which is matched */
//                                                                        /* with y, to the queue */
//                    add_to_tree(yx[y], x);                              /* add edges (x,y) and (y,yx[y]) to the tree */
//                }
//            if (y < n) break;                                           /* augmenting path found! */
//        }
//        if (y < n) break;                                               /* augmenting path found! */
//
//         /* augmenting path not found, so improve labeling. "update_labels" */
//		for (j = 0; j < n; j++)                                         /* calculate delta using slack */
//        if (!Th[j])
//            delta = whichmin(delta, slack[j]);
//        for (i = 0; i < n; i++)                                         /* update X labels */
//            if (Sh[i]) lx[i] -= delta;
//        for (j = 0; j < n; j++)                                         /* update Y labels */
//            if (Th[j]) ly[j] += delta; 
//        for (j = 0; j < n; j++)                                         /* update slack array */
//            if (!Th[j])
//                slack[j] -= delta;
//
//        wr = rd = 0;                
//        for (y = 0; y < n; y++)        
//        /* in this cycle we add edges that were added to the equality graph as a
//        result of improving the labeling, we add edge (slackx[y], y) to the tree if
//        and only if !Th[y] &&  slack[y] == 0, also with this edge we add another one
//        (y, yx[y]) or augment the matching, if y was exposed */
//            if (!Th[y] &&  slack[y] == 0)
//            {
//                if (yx[y] == -1)                                        /* exposed vertex in Y found - augmenting path exists! */
//                {
//                    x = slackx[y];
//                    break;
//                }
//                else
//                {
//                    Th[y] = 1;                                          /* else just add y to Th, */
//                    if (!Sh[yx[y]])    
//                    {
//                        q[wr++] = yx[y];                                /* add vertex yx[y], which is matched with */
//                                                                        /* y, to the queue */
//                        add_to_tree(yx[y], slackx[y]);                  /* and add edges (x,y) and (y, */
//                                                                        /* yx[y]) to the tree */
//                    }
//                }
//            }
//        if (y < n) break;                                               /* augmenting path found! */
//    }
//
//    if (y < n)                                                          /* we found augmenting path! */
//    {
//        max_match++;                                                    /* increment matching */
//        /* in this cycle we inverse edges along augmenting path */
//
//        for (cx = x, cy = y; cx != -2; cx = prev[cx], cy = ty)
//        {
//            ty = xy[cx];
//            yx[cy] = cx;
//            xy[cx] = cy;
//        }
//
//        augment();                                                      /* recall function, go to step 1 of the algorithm */
//    }
//}
//

/*----------------------------------------------------------------------------------------------------*/
int hungarian(void);
//{
//    int ret = 0;                      /* weight of the optimal matching */
//    int x, y;
//
//    max_match = 0;                    /* number of vertices in current matching */
//    memset(xy, -1, sizeof(xy));    
//    memset(yx, -1, sizeof(yx));
//
//    memset(lx, 0, sizeof(lx));        /* step 0, "init_labels" */
//    memset(ly, 0, sizeof(ly));
//
//    for (x = 0; x < n; x++)
//        for (y = 0; y < n; y++)
//            lx[x] = whichmax(lx[x], cost[x][y]);                    
//
//    augment();                        /* steps 1-3 */
//    for (x = 0; x < n; x++)       /* forming answer there */
//        ret += cost[x][xy[x]];
//
//    return (ret);
//}


/*----Calculate the partition distance between A and B which has Ne elements--------------------------*/
int Partdis ( const int Ne, const int *A, const int *B );
//{
//	int ii, jj, kk, Pos, sum, Dis;
//	int *Abs=NULL, *Bbsbar=NULL;
//
//	/* n represents the number of clusters in A and B */
//	for (ii=0, n=0; ii<Ne; ii++)
//	{
//		if ( n < (int) A[ii] )
//			n = (int) A[ii];
//		if ( n < (int) B[ii] )
//			n = (int) B[ii];
//	}
//	n += 1; /* because the cluster number in A and B starts with 0 */
//
//	if ( n > Nv )
//	{
//		printf ( "Number of clusters in sampled partitions is larger than %d\n", Nv);
//		return (0);
//	}
//
//	Abs = calloc ( Ne * n, sizeof ( int ) );    if(Abs==NULL) alert();
//	Bbsbar = calloc ( Ne * n, sizeof ( int ) ); if(Bbsbar==NULL) alert();
//
//	for (ii=0; ii<n; ii++)
//		for (jj=0; jj<Ne; jj++ )
//		{
//			Pos = ii * Ne + jj;
//			if ( A[jj] == (int) ii ) Abs [ Pos ] = 1;
//			if ( B[jj] != (int) ii ) Bbsbar [ Pos ] = 1;
//		}
//
//	for (ii=0; ii<n; ii++ )
//	{
//		for (jj=0; jj<n; jj++ )
//		{
//			sum=0;
//			for ( kk=0; kk<Ne; kk++ )
//				sum += ( Abs [ jj*Ne+kk ] * Bbsbar [ ii*Ne+kk ] );
//
//			cost [ii][jj] = -sum;
//		}
//	}
//
//	Dis = hungarian ();
//	Dis *= -1;
//
//	free(Abs); Abs=NULL; free(Bbsbar);  Bbsbar=NULL;
//
//	return(Dis);
//
//}

/*----Main ( calculate the mean partition distance )--------------------------------------------------*/
void Meanpartdis ( int *sampledz, int ne, int ns, int *meanpartition);
//{
//	int ii, jj;
//	int *Meanp, *propMeanp, kmeanp=0, kk, min;
//	int Distsum, propDistsum=0, ele=0, Move=0, Notmove=0;
//	int poplabel, labelupdate;
//
//	// VS
//	// ns is the number of sampled assignment vectors
//	// Note that in this case the assignment vector is treated has
//	// a single array, instead of being seen as a 2D array (i.e. a matrix)
//	Distsum = ne * ns;
//
//	// ne is the number of elements
//	// meanp has size ne
//	Meanp = malloc ( sizeof ( int ) * ne); if(Meanp==NULL) alert();
//	// propMeanp also has size ne
//	propMeanp = malloc ( sizeof ( int ) * ne); if(propMeanp==NULL) alert();
//	
//	// is this initializing the meanp with the firt sampledz
//	for (ii=0; ii<ne; ii++ ){ Meanp[ii] = sampledz [ii];}
//
//	
//	do
//	{
//		kmeanp = 0;
//		for (ii=0; ii<ne; ii++ )
//		{
//			if ( Meanp[ii] > kmeanp )
//				kmeanp = Meanp [ii];
//			propMeanp [ii] = Meanp [ii];
//		}
//		kmeanp += 1;	
//
//		for (kk=0; kk<=kmeanp; kk++ )
//		{
//			propMeanp[ele] = kk;
//			
//			propDistsum = 0;
//			for (jj=0; jj<ns; jj++ )
//			{
//				propDistsum += Partdis (ne, propMeanp, sampledz + ( jj * ne ) );
//
//			}
//			
//			if ( propDistsum < Distsum )
//			{
//				Distsum = propDistsum;
//				for (jj=0; jj<ne; jj++ ) { Meanp[jj] = propMeanp[jj]; }
//				Move += 1;
//				Notmove = 0;
//			}
//		}
//
//		ele += 1;
//
//		if ( Move == 0 ) { Notmove += 1;} else { Move = 0; }
//		if ( ele == ne )
//			ele = 0;
//
//	} while ( Notmove<( ne - 1) );
//
//
//	min = Meanp [0];
//	for ( ii=1; ii<ne; ii++ )
//		if ( Meanp [ii] < min )
//			min = Meanp [ii];
//
//	jj=0; poplabel=1;
//	do
//	{
//		labelupdate = 0;
//		for ( ii=0; ii<ne; ii++ )
//			if ( Meanp [ii] == min )
//			{
//				meanpartition [ii] = poplabel;
//				jj ++;
//				labelupdate ++;
//			}
//
//		min ++;
//
//		if (labelupdate>0)
//			poplabel ++;
//
//	} while ( jj<ne );
//
//
//	
//	free(propMeanp); free(Meanp);
//
//}


/* Calculate the co-assignment probs ---------------------------------------------------------------------------------------------------*/
void Coassignprobs ( int *sampledz, int ni, int ns, double *coaprob );
//{
//	int ii, jj, ss, rr;
//
//    for ( ss=0; ss<ns; ss++ )
//	{
//		rr = 0;
//		for ( ii=0; ii<(ni-1); ii++ )
//		{
//			for ( jj=ii+1; jj<ni; jj++ )
//			{
//				if ( sampledz [ ss * ni + ii ] == sampledz [ ss * ni + jj ] )
//					coaprob [rr] += 1.0;
//				rr ++;
//			}
//		}
//	}
//
//	for ( ii=0; ii<rr; ii++ )
//		coaprob [ii] /= (double) ns;
//
//}
//


/*----alert when allocation of memory is failed-----------------------------------------------------------------------------------*/
void alert ( void );
/*{
	puts ( "memory allocation failure. try again" );
}*/



// INDPROB
// returns the uncertainty probability of each individual being assigned
// to the cluster identified in the meanpartition.
// Function coded by VS
void indprob (int *meanpart, int *sampledz, int ni, int ns, double *returnprob);

// FACTORIAL
// returns the factorial of an integer smaller than 13
int factorial (int n);
