#include<iomanip>
#include <mpi.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <string.h>
#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <time.h>

extern "C++" void sort( long *A, long q, long r );

using namespace std;

/*  globals */
int numnodes,rank,mpi_err, cilk_err;
#define mpi_root 0
/* end globals  */



int binaryIndexSearch(long *, long, long , long );

void init_it(int  *argc, char ***argv);

void init_it(int  *argc, char ***argv) {
	mpi_err = MPI_Init(argc,argv);
    mpi_err = MPI_Comm_size( MPI_COMM_WORLD, &numnodes );
    mpi_err = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
}



int main(int argc,char *argv[]){
	long *part_array,*send_array,*recv_array, *output_array;
	long *sub_sorted, *sort_keys, *pivot_elements;
	long *global_pivots, *index_pivots, *index_to_send, *length_to_send, *length_to_receive;
	long *recv_buffer, *prefix_length, *prefix_length_to_send, *index_to_receive;
	long *length_from_all_processor, *prefix_final;
	int *send_length_array, *displs;
	long count, temp;
	long size,mysize,i,k,j,total,p;
	long q, partition;
	long recv_length=0;
	long add_to_size=0;
	long remainder=0;
	string line;
	long *input_array;	
	
/*	time_t t0, t1;
        clock_t c0,c1;
*/

	MPI_Status status, stat[6]; 					// status of MPI_Recv
	MPI_Request request, send_request[12], recv_request[12]; 				// capture request of MPI_Isend
	
	ifstream myfile(argv[1]);
	if(myfile.is_open())
		{
			if(myfile.good())
			{
				for(i=0; !myfile.eof();i++)
				{
					if(i==0)
					{
						getline(myfile, line);
		  				size = atol(line.c_str());
						input_array = new long[size];
					}
					else
					{
						getline(myfile, line);
						input_array[i-1]=atol(line.c_str());
					}
				}
			}

		}

	k=5;
	init_it(&argc,&argv);
	q= k*numnodes;
	partition = q-1;
	
	p = numnodes-1;
	
/* each processor will get count elements from the root */
	
	
	if(rank == mpi_root)
	{
		send_array = input_array;
	}	

	mpi_err = MPI_Barrier( MPI_COMM_WORLD);
	count=size/numnodes;
	remainder = size%numnodes;
	
	if(rank == (numnodes-1))
	{
		count = count+remainder;
	}
	
	
	part_array=(long*)malloc(count*sizeof(long));
	
	mpi_err = MPI_Barrier( MPI_COMM_WORLD);

/* create the data to be sent on the root */
	if(rank == mpi_root){

//		printf("%ld \n",(count/q));
//		printf("Size is = %ld\n",size);
		recv_array=(long *)malloc(size*sizeof(long));
		output_array=(long *)malloc(size*sizeof(long));
		memset(output_array,0,size*sizeof(long));
/*		for(i=0;i<size;i++)
	    		printf("%ld ",send_array[i]);
		printf("\n The input output line \n");
*/
		
		sort_keys = (long *)malloc(numnodes * partition * sizeof(long));
		
		send_length_array = (int *)malloc(numnodes*sizeof(int));
		for(i=0;i<numnodes;i++)
		{
			if(i != (numnodes-1))
				send_length_array[i]=(int)count;
			else
				send_length_array[i]=(int)(count+remainder);
		}
		
		displs = (int *)malloc(numnodes*sizeof(int));
		for(i=0;i<numnodes;i++)
		{
			displs[i]=(int)(i*count);
		}
		
	}
	
/*	if(rank == mpi_root)
	{
		for(i=0;i<numnodes;i++)
	    		printf("%d ",send_length_array[i]);
		printf("\n The length info \n");
		for(i=0;i<numnodes;i++)
	    		printf("%d ",displs[i]);
		printf("\n The displacement \n");
	}
*/	
	
	mpi_err = MPI_Barrier( MPI_COMM_WORLD);
	
/* send different data to each processor */
/*	mpi_err = MPI_Scatter(	send_array, count,   MPI_LONG,
						    part_array,    count,   MPI_LONG,
	                 	    mpi_root,
	                 	    MPI_COMM_WORLD);
*/
	mpi_err = MPI_Scatterv(send_array, send_length_array, displs, MPI_LONG,
							part_array, count, MPI_LONG,
							mpi_root,
							MPI_COMM_WORLD);
	

/*	t0=time(NULL);
        c0=clock();
	 printf ("\tbegin (wall):            %ld\n", (long) t0);
        printf ("\tbegin (CPU):             %d\n", (int) c0);
*/	

/* each processor does a local sum */
	total=0;



//start your code here shared memory quick sort
/*	for(i=0;i<count;i++)
	{
	    temp = part_array[i];
		j = i-1;
		while(temp<part_array[j] && j>=0)
		{
			part_array[j+1] = part_array[j];
			j = j-1;
		}
		part_array[j+1] = temp;
	}
*/

	sort(part_array, 0, count-1);

	
//end your code here
	
	if( (count/q) >= 1)
	{
		
		sub_sorted = (long *)malloc(partition*sizeof(long));
		for(i=0; i<partition; i++)
		{
			sub_sorted[i] = part_array[(i+1)*(count/q)];
		}
		
	}
	mpi_err = MPI_Gather(sub_sorted, partition,  MPI_LONG, 
						sort_keys,partition,  MPI_LONG, 
	                 	mpi_root,                  
	                 	MPI_COMM_WORLD);

						
	if( rank == mpi_root)
	{
	
		// start your code here (shared memory quick sort 
/*		for(i=0;i< (numnodes*partition);i++)
		{
			temp = sort_keys[i];
			j = i-1;
			while(temp<sort_keys[j] && j>=0)
			{
				sort_keys[j+1] = sort_keys[j];
				j = j-1;
			}
			sort_keys[j+1] = temp;
		}
*/
		sort(sort_keys, 0, (numnodes*partition)-1);


		// end you code here
/*		for(i=0;i< (numnodes*partition);i++)
			printf("%ld ",sort_keys[i]);
		printf("\n Sorted partial keys \n ");
*/		
		p= numnodes-1;
		if((partition) >= 1)
		{
			pivot_elements = (long *)malloc(p*sizeof(long));
			for(i=0; i<p; i++)
			{
				pivot_elements[i] = sort_keys[(i+1)*partition];
			}
		}
/*		for(i=0;i<p;i++)
			printf("%ld ",pivot_elements[i]);
		printf("\n Pivot elements \n");
*/
	}

	mpi_err = MPI_Barrier( MPI_COMM_WORLD);
	//Code for broadcasting the pivot_elements
	p = numnodes-1;
	global_pivots = (long *)malloc(p*sizeof(long));
	memset(global_pivots, 0, p*sizeof(long));
	if(rank == mpi_root)
	{
		global_pivots = pivot_elements;
	}

	mpi_err = MPI_Barrier( MPI_COMM_WORLD);

	mpi_err = MPI_Bcast(global_pivots, p, MPI_LONG, mpi_root, MPI_COMM_WORLD);
/*	if(rank == (mpi_root+1))
	{
		for(i=0;i<p;i++)
			printf("%d ",global_pivots[i]);
		printf("\n Global Pivot elements \n");
	}
*/
	
	index_pivots = (long *)malloc(p*sizeof(long));
	memset(index_pivots, 0, p*sizeof(long));
	//Getting indexes of the global pivots in the local sorted lists
	for(i=0; i<p; i++)
	{
		index_pivots[i] = binaryIndexSearch(part_array, global_pivots[i], 0, count-1);
	}
	
	
	
	//Now will set the lengths for the proccessor i to get the ith bucket
	length_to_send = (long *)malloc(numnodes*sizeof(long));
	memset(length_to_send, 0, numnodes*sizeof(long));
	for(i =0; i<numnodes; i++)
	{
		if(i == 0)
		{
			length_to_send[i] = index_pivots[i]+1;
		}
		if( i > 0 && i <(numnodes-1))
		{
			if(index_pivots[i-1] == index_pivots[i])
			{
				length_to_send[i]=0;
			}
			else
			{
				length_to_send[i]=index_pivots[i] - index_pivots[i-1];
			}
		}
		if( i == (numnodes-1))
		{
			length_to_send[i]=(count-1) - index_pivots[i-1];
		}
	}
	
	
	prefix_length_to_send = (long *)malloc(numnodes*sizeof(long));
	memset(prefix_length_to_send, 0, numnodes*sizeof(long));
	for(i=0; i<numnodes;i++)
	{
		if(i==0)
			prefix_length_to_send[i]=length_to_send[i];
		else
			prefix_length_to_send[i]=prefix_length_to_send[i-1]+length_to_send[i];
	}
	for(i=0; i<numnodes;i++)
	{
		prefix_length_to_send[i] = prefix_length_to_send[i]-1;
	}
	index_to_send = (long *)malloc(numnodes*sizeof(long));
	memset(index_to_send, 0, numnodes*sizeof(long));
	for(i=0;i<numnodes;i++)
	{
		if(i==0)
			index_to_send[i]=0;
		else
		{
			if(prefix_length_to_send[i] != prefix_length_to_send[i-1])
				index_to_send[i]= prefix_length_to_send[i-1]+1;
		}
	}
	
	mpi_err = MPI_Barrier( MPI_COMM_WORLD);

	//Now will send the length to each nodes
	length_to_receive = (long *)malloc(numnodes*sizeof(long));
	memset(length_to_receive, 0, numnodes*sizeof(long));
//	for(i=0; i< numnodes; i++)
/*	{
		
		mpi_err = MPI_Isend(&length_to_send[0], 1, MPI_LONG, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &send_request[0]);
		mpi_err = MPI_Irecv(&length_to_receive[0], 1, MPI_LONG, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &recv_request[0]);
	}
	
//	for(i=0; i< numnodes; i++)
	{
		mpi_err = MPI_Wait(&recv_request[0], &status);
		mpi_err = MPI_Wait(&send_request[0], &status);
	}
*/

	for(i=0; i<numnodes;i++)
	{
		mpi_err = MPI_Scatter(length_to_send, 1,   MPI_LONG,
								&length_to_receive[i], 1,   MPI_LONG,
								i,
								MPI_COMM_WORLD);
							
	}						
	mpi_err = MPI_Barrier( MPI_COMM_WORLD);
	for(i=0; i<numnodes;i++)
	{
		recv_length = recv_length + length_to_receive[i];
	}
	
	
	//Now will send that portion to each processor.
	recv_buffer = (long *)malloc(recv_length*sizeof(long));
	memset(recv_buffer,0,recv_length*sizeof(long));
	
	prefix_length = (long *)malloc(numnodes*sizeof(long));
	for(i=0; i<numnodes;i++)
	{
		if(i==0)
			prefix_length[i]=length_to_receive[i];
		else
			prefix_length[i]=prefix_length[i-1]+length_to_receive[i];
	}
	for(i=0; i<numnodes;i++)
	{
		prefix_length[i] = prefix_length[i]-1;
	}
	index_to_receive = (long *)malloc(numnodes*sizeof(long));
	memset(index_to_receive, 0, numnodes*sizeof(long));
	for(i=0;i<numnodes;i++)
	{
		if(i==0)
			index_to_receive[i]=0;
		else
		{
			index_to_receive[i]= prefix_length[i-1]+1;
		}
	}
	
/*	if(rank == 1)
		mpi_err = MPI_Isend(part_array, length_to_send[0], MPI_LONG, 0, 2, MPI_COMM_WORLD, &send_request[0]);
	if(rank == 0)
		mpi_err = MPI_Irecv(&recv_buffer[length_to_receive[0]], length_to_receive[1], MPI_LONG, 1, 2, MPI_COMM_WORLD, &recv_request[0]);
*/
/*	if(rank == 1)
		mpi_err = MPI_Send(part_array, length_to_send[0], MPI_LONG, 0, 2, MPI_COMM_WORLD);
	
	if(rank == 4)
	{
		mpi_err = MPI_Send(part_array, length_to_send[0], MPI_LONG, 0, 4, MPI_COMM_WORLD);
	}
	if(rank == 5)
	{
		mpi_err = MPI_Send(part_array, length_to_send[0], MPI_LONG, 0, 5, MPI_COMM_WORLD);
	}
	if(rank == 6)
	{
		mpi_err = MPI_Send(part_array, length_to_send[0], MPI_LONG, 0, 6, MPI_COMM_WORLD);
	}
	if(rank == 7)
	{
		mpi_err = MPI_Send(part_array, length_to_send[0], MPI_LONG, 0, 7, MPI_COMM_WORLD);
	}
	if(rank == 8)
	{
		mpi_err = MPI_Send(part_array, length_to_send[0], MPI_LONG, 0, 8, MPI_COMM_WORLD);
	}
	if(rank == 9)
	{
		mpi_err = MPI_Send(part_array, length_to_send[0], MPI_LONG, 0, 9, MPI_COMM_WORLD);
	}
	if(rank == 10)
	{
		mpi_err = MPI_Send(part_array, length_to_send[0], MPI_LONG, 0, 10, MPI_COMM_WORLD);
	}
	if(rank == 0)
	{
		mpi_err = MPI_Recv(&recv_buffer[prefix_length[0]+1], length_to_receive[1], MPI_LONG, 1, 2, MPI_COMM_WORLD, &status);
		mpi_err = MPI_Recv(&recv_buffer[prefix_length[3]+1], length_to_receive[4], MPI_LONG, 4, 4, MPI_COMM_WORLD, &status);
		mpi_err = MPI_Recv(&recv_buffer[prefix_length[4]+1], length_to_receive[5], MPI_LONG, 5, 5, MPI_COMM_WORLD, &status);
		mpi_err = MPI_Recv(&recv_buffer[prefix_length[5]+1], length_to_receive[6], MPI_LONG, 6, 6, MPI_COMM_WORLD, &status);
		mpi_err = MPI_Recv(&recv_buffer[prefix_length[6]+1], length_to_receive[7], MPI_LONG, 7, 7, MPI_COMM_WORLD, &status);
		mpi_err = MPI_Recv(&recv_buffer[prefix_length[7]+1], length_to_receive[8], MPI_LONG, 8, 8, MPI_COMM_WORLD, &status);
		mpi_err = MPI_Recv(&recv_buffer[prefix_length[8]+1], length_to_receive[9], MPI_LONG, 9, 9, MPI_COMM_WORLD, &status);
		mpi_err = MPI_Recv(&recv_buffer[prefix_length[9]+1], length_to_receive[10], MPI_LONG, 10, 10, MPI_COMM_WORLD, &status);
	}
	*/
	
	mpi_err = MPI_Barrier( MPI_COMM_WORLD);
	for(j=0; j<numnodes;j++)
	{
		for(i =0; i<numnodes; i++)
		{
			if(j!=i)
			{
				if(rank == i)
				{
//					printf("\nSending data by rank = %d\n",i);
					mpi_err = MPI_Send(&part_array[index_to_send[j]], length_to_send[j], MPI_LONG, j, i, MPI_COMM_WORLD);
//					mpi_err = MPI_Isend(part_array, length_to_send[j], MPI_LONG, j, i, MPI_COMM_WORLD, &request);
					
				}
				if(rank == j)
				{
//					printf("\n Waiting to recieve data from rank=%d \n",i);
					mpi_err = MPI_Recv(&recv_buffer[index_to_receive[i]], length_to_receive[i], MPI_LONG, i, i, MPI_COMM_WORLD, &status);
//					mpi_err = MPI_Irecv(&recv_buffer[prefix_length[i-1]+1], length_to_receive[i], MPI_LONG, i, i, MPI_COMM_WORLD, &request);
				}
			}
		}
	}
	mpi_err = MPI_Barrier( MPI_COMM_WORLD);
	if((length_to_send[rank] == length_to_receive[rank]) && (length_to_send[rank]>0))
	{
		memcpy(&recv_buffer[index_to_receive[rank]], &part_array[index_to_send[rank]], length_to_receive[rank]*sizeof(long));
	}

/*	if((length_to_send[rank] == length_to_receive[rank]) && (length_to_send[rank]>0))
	{
		for(i=0;i<length_to_receive[rank];i++)
		{
			memcpy(&recv_buffer[index_to_receive[rank]+i], &part_array[index_to_send[rank]+1], sizeof(int));
		}
	}
*/
//	mpi_err = MPI_Wait(&recv_request[0],&status);
//	mpi_err = MPI_Wait(&send_request[0],&status);
	
	mpi_err = MPI_Barrier( MPI_COMM_WORLD);
	//Now will sort the recv array
	// start your code here (shared memory quick sort 
/*		for(i=0;i<recv_length;i++)
		{
			temp = recv_buffer[i];
			j = i-1;
			while(temp<recv_buffer[j] && j>=0)
			{
				recv_buffer[j+1] = recv_buffer[j];
				j = j-1;
			}
			recv_buffer[j+1] = temp;
		}
*/
	sort(recv_buffer, 0, (recv_length)-1);


	// end you code here

	mpi_err = MPI_Barrier( MPI_COMM_WORLD);
	length_from_all_processor = (long *)malloc(numnodes*sizeof(long));
	memset(length_from_all_processor, 0, numnodes*sizeof(long));
	mpi_err = MPI_Gather( &recv_length, 1, MPI_LONG,
							length_from_all_processor, 1, MPI_LONG,
							mpi_root,
							MPI_COMM_WORLD);
	


	prefix_final = (long *)malloc(numnodes*sizeof(long));
	memset(prefix_final, 0, numnodes*sizeof(long));
	for(i=0; i<numnodes;i++)
	{
		if(i==0)
			prefix_final[i]=length_from_all_processor[i];
		else
			prefix_final[i]=prefix_final[i-1]+length_from_all_processor[i];
	}
	
	mpi_err = MPI_Barrier( MPI_COMM_WORLD);
	
/*	
	t1=time(NULL);
        c1=clock();
        printf ("\telapsed wall clock time: %ld\n", (long) (t1 - t0));
        printf ("\telapsed CPU time:        %f\n", (float) (c1 - c0)/CLOCKS_PER_SEC);
*/

	for(i =1; i<numnodes; i++)
	{
		
		if(rank == i)
		{
			mpi_err = MPI_Send(recv_buffer, recv_length, MPI_LONG, 0, i, MPI_COMM_WORLD);
			
		}
		if(rank == 0)
		{
			mpi_err = MPI_Recv(&output_array[prefix_final[i-1]], length_from_all_processor[i], MPI_LONG, i, i, MPI_COMM_WORLD, &status);

		}
		
	}
	
	mpi_err = MPI_Barrier( MPI_COMM_WORLD);

	if(rank == mpi_root)
	{
		memcpy(output_array, recv_buffer, recv_length*sizeof(long));
	}
	
	mpi_err = MPI_Barrier( MPI_COMM_WORLD);

 /*   mpi_err = MPI_Gather(part_array,   count,  MPI_LONG, 
						recv_array,count,  MPI_LONG, 
	                 	mpi_root,                  
	                 	MPI_COMM_WORLD);
*/
/* the root prints the global sum */
	if(rank == mpi_root){

/*		for(i=0;i<size;i++)
			printf("%ld ",recv_array[i]);
		printf("\n results from all processors\n");
		for(i=0;i<p;i++)
			printf("%ld ",index_pivots[i]);
		printf("\n Index of the sorted part arrays acc to global pivot\n \n");
		for(i=0;i<numnodes;i++)
			printf("%ld ",length_to_send[i]);
		printf("\n Length details\n \n");
		for(i=0;i<numnodes;i++)
			printf("%ld ",length_to_receive[i]);
		printf("\n Length details received\n \n");
		for(i=0;i<numnodes;i++)
			printf("%ld ",prefix_length[i]);
		printf("\n Prefix Length details received\n \n");
		
		for(i=0; i<recv_length;i++)
			printf("%ld ",recv_buffer[i]);
		printf("\n Receive Array\n\n");
		for(i=0; i<recv_length;i++)
			printf("%ld ",length_from_all_processor[i]);
		printf("\n Receive Final Array Length\n\n");
		for(i=0;i<numnodes;i++)
			printf("%ld ",prefix_final[i]);
		printf("\n Prefix Length details received\n \n");
*/		for(i=0;i<size;i++)
			printf("%ld \n",output_array[i]);
//		printf("\n results from all processors\n\n");
		
	}

/*	if(rank == 11)
	{
		for(i=0;i<p;i++)
			printf("%d ",index_pivots[i]);
		printf("\n Index of the sorted part arrays acc to global pivot\n \n");
		for(i=0;i<numnodes;i++)
			printf("%d ",length_to_send[i]);
		printf("\n Length details\n \n");
		for(i=0;i<numnodes;i++)
			printf("%d ",prefix_length_to_send[i]);
		printf("\n Prefix Length details sent\n \n");

		for(i=0;i<numnodes;i++)
			printf("%d ",length_to_receive[i]);
		printf("\n Length details received\n \n");
		for(i=0;i<numnodes;i++)
			printf("%d ",prefix_length[i]);
		printf("\n Prefix Length details received\n \n");
		for(i=0;i<recv_length;i++)
			printf("%d ",recv_buffer[i]);
		printf("\n Receive Array for rank= %d\n \n", rank);
		printf("\n Receive length is : %d\n \n",recv_length);
		printf("\n length to receive= %d",length_to_receive[rank]);
	}
*/	
//	printf("Rank is %d  and length to receive is %d \n", rank, recv_length);

	mpi_err = MPI_Barrier( MPI_COMM_WORLD);
    mpi_err = MPI_Finalize();
}

int binaryIndexSearch(long *A, long key, long low, long high)
{
  if (high < low)
    return (low-1);
  else
    {
      int mid=(low+high)/2;
      if (A[mid] > key)
        return binaryIndexSearch(A, key, low, mid-1);
      else if (key >= A[mid])
        return binaryIndexSearch(A, key, mid+1, high);
      else
        return mid;
    }
}
