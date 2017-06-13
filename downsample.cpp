
/** 
 *  @file        downsample.cpp
 *  @author  Chandrahas J R
 *  @date      5/21/2017  
 *  @version 1.0 
 *  
 *  @brief Programming assignment for Seung Lab 
 *
 *  @section DESCRIPTION
 *  
 *  This is a program that let's a user create a 2D array of size dimA X dimB where dimA and dimB are powers of  2.
 *  
 * The program then prints all the downsampled verison's of the 2D array where the dimensions of the array reduce by 2^l where l = {1....min(dimA, dimb)} 
 *
 */


#include <iostream>
#include <boost/multi_array.hpp>
#include <cassert>
#include <math.h>
#include <vector>
#include <map>
#include <cstdlib>
#include <thread>
#include <future>
#include <list>
#include <unordered_map>
#include <mutex>
#include <chrono>


using namespace std;
using namespace std::chrono;

//Set the max threads that the program can create, to perform downsampling
//Best value is the number of cores in the machine 
#define MAX_THREADS 4

typedef boost::multi_array< unsigned int, 2> baseImage;
typedef baseImage::index index;
baseImage::extent_gen extents;


std::mutex mapMutexThread1;
std::mutex mapMutexThread2;
std::mutex mapMutexThread3;
std::mutex mapMutexThread4;

/*----------------------------------------------------------
* DESCRIPTION
* 
* This class manages the information about a 2x2 block of data. 
* 
* member variables m_mode and m_count store the mode value and its count .
*
* For example, for the 2x2 block shown below,
*[  1,  1]
*[  1,  2]
* m_mode = 1
* m_count = 3
*
* m_depth indicates how many recursive calls to startThreading happened befoe the object was createad.
* ---------------------------------------------------------------
*/
class modeMap
{
	private:
		std::map< unsigned int, unsigned int> cube;
		unsigned int m_mode;
		unsigned int m_count;
		int m_threadNumber;
		int m_depth;

     
	public:
    /**
        Constructor
        */
		modeMap(): m_mode(0), m_count(0)  
		{
		}
		
        size_t getMapSize() { return cube.size(); }
		unsigned int getMode() { return m_mode; }

		unsigned int getCount(){ return m_count; }

         /**
        
        Add's a new element to the map if it isn't present already.
        Increases the value by one to count for the element in the 2x2 block.
        
        @param data - element to be added the map
        @return void      
       */
		void addElement(unsigned int data, int count = 1)
		{	
			cube[data]  += count;
            return;
		}

        
        /**
        
        Function to calcuate the mode and count of the 2x2 block. 
        This function is called after the cube data member has been populated 
        
        @param void
        @return void   
        
        */
        
		void calculateMode()
		{		
			for(auto element = cube.begin(); element != cube.end(); ++element)
			{	
				if(element->second > m_count) 
				{	
					m_mode = element->first;
					m_count = element->second;
				}
			}
		}

        
        // Setter function to set the depth of the cube in the base image
		void setDepth(int depth) { m_depth = depth; }
	
        // Setter function to set the threadNumber by which the cube was read by
		void setthreadNumber(int threadNumber) { m_threadNumber = threadNumber; }

        //Overloaded operator used to merge modeMap objects while downsampling.
        friend modeMap operator+ (modeMap input1, modeMap input2)
        {
            map< unsigned int, unsigned int> inMap1 = input1.getCube();
            
            for(auto iter = inMap1.begin(); iter != inMap1.end(); ++iter)
            {
                input2.addElement(iter->first,iter->second);
            }
            
            return input2;
        
        }
        
        //Getter function required for downsampling. 
        map< unsigned int, unsigned int> & getCube()
        {
            return cube;
        }
        
		~modeMap()
		{
		}

}; 


/*----------------------------------------------------------
* DESCRIPTION
* 
* This class manages the information for an entire 2D Array of data. 
* 
* Member variables are described as below.
*
* dimA - Number of Rows
* dimB - Number of Cols
* minDim - min (dimA, dimB)
* m_numCols - number of Cols of 2X2 blocks 
* m_depth -  the number of recursive calls to startThreading happened befoe the object was createad.
* m_mapCollectThread1,m_mapCollectThread1,m_mapCollectThread1,m_mapCollectThread4 - A map of modeMap objects in a sub matrix of the baseImage. From startThreading function, a minimum of 4 threads will be created and each,
                                                                                                                                          one these maps store the mode information for each small 2x2 block within that submatrix. The key for each modeMap in these maps is the location of the 2X2 blocks in the entire baseImage.
* m_globalMap - Once all the above maps hav been created, they are stored together in one map for the enitre baseImage
* m_baseImage is a 2 dimensional boost Multi Array

For Example,

Given a baseImage of size, 8x8 the image is divided into 4 threads and each threads creates 4 more untill the total number of threads created reaches the max number defined above.

3 4 0 6 | 6 4 3 1
1 4 4 4 | 0 1 4 0
6 1 4 3 | 7 0 4 3
4 5 2 2 | 1 4 7 4
--------------------
7 5 8 4 | 1 2 4 0
4 8 4 4 | 1 7 3 7
8 7 8 4 | 7 3 8 2
6 8 2 7 | 1 1 1 8

m_mapCollectThread1(2,3,4) save data about each 2x2 block in the baseImage and is grouped based on which thread read that 2x2 block.
If max threads = 4, then each m_mapCollectThread1(2,3,4) would have 4 maps each with m_mode, m_count.
The m_globalMap would later have all the 16 maps corresponding to the 16 2x2 blocks in the base image. The index for each 2X2 block is calculated as named in the following convention 

01 02 | 03 04
05 06 | 07 08
-----------------
09 10 | 11 12
13 14 | 15 16

Once we have the data for a 2x2 block, we can further derive the downsampled image by grouping these 2x2 blocks together. 
 
* ---------------------------------------------------------------
*/

class twoDArray
{
        private:
            int m_dimA, m_dimB, m_minDim, m_depth;
            
            int m_downRows, m_downCols;
            
            map < int, modeMap > m_mapCollectThread1;
            map < int, modeMap > m_mapCollectThread2;
            map < int, modeMap > m_mapCollectThread3;
            map < int, modeMap > m_mapCollectThread4;
            map < int, modeMap > m_globalMap;
            
            baseImage m_baseImage;
            
            public:
            
            /**
            
            Constructor for the class. Initializes all member variables and also fills the matrix with random values from 1-8.           
              
            */
            twoDArray(int dimA, int dimB) : m_dimA(dimA), m_dimB(dimB), m_minDim((m_dimA > m_dimB) ? m_dimB : m_dimA), m_depth(0), m_baseImage(boost::extents[m_dimA][m_dimB]), m_downRows(dimA/2), m_downCols(dimB/2)

            {
                while(m_minDim != 1)
                {	
                    m_minDim /= 2;
                    ++m_depth;
                } 

                for(index i = 0; i < m_dimA; ++i)
                {
                    for(index j = 0; j < m_dimB; ++j)
                    {
                        m_baseImage[i][j] = rand()%9;
                        
                        //uncommment the below TWO statements to see the original matrix
                        //std::cout << m_baseImage[i][j] << " ";
                    }
                    //std::cout << std::endl;
                }
                
                //Restore min Dimension. We will need this when we are printing downsampled versions
                m_minDim  = (m_dimA > m_dimB) ? m_dimB : m_dimA;
                //m_minDim /= 2;
            }
            
            /**
            
            Function to create a modeMap object for every 2x2 block. The function add's all the four elements first and then calls the caluclateMode(), setDepth and setThreadNumber function.             
             
             @param 
                    rowStart, rowEnd - row positions of the elements in the baseImage 
                    colStart, colEnd - col psitiiop pf tje elemtnes in the baseImage.
                    depth - Depth of the cube in the baseImage 
                    threadNumber - thread reading a particular 2x2 block in the baseImage.

             @return 
                    modeMap object with all the elements, mode and count data.            
            */
            
            
            modeMap findMode(int rowStart, int rowEnd, int colStart, int colEnd, int depth, int threadNumber)
            {
 
                modeMap result;
               
                for(int i = rowStart; i < rowEnd; ++i)
                {	
                    for(int j = colStart; j < colEnd; ++j)
                    {   
                        unsigned int currentElement = m_baseImage[i][j];
                        result.addElement(currentElement);
                    }
                }
                
                result.setDepth(depth);
                result.setthreadNumber(threadNumber);	
                result.calculateMode();
                return result;

            }            
            
             /**
            
            Function to add a mapMode object to one of the m_mapCollectThread1(2,3,4) objects.            
            Mutexes are added to prevent race aroud conditions
            
             @param 
                    result - modeMap object with all the elements, mode and count data.
                    index - index of the 2x2 block in the baseImage.
                    threadNumber - thread reading a particular 2x2 block in the baseImage.

             @return 
                    void           
            */
            
            
            void addBlock(int threadNumber, int index, modeMap result) 
            {
               // std::cout << "I'm Inside add Block " << std::endl;
                switch(threadNumber){

                    case 1 :
                        { 
                        std::lock_guard<std::mutex> guard(mapMutexThread1);
                        m_mapCollectThread1.insert(std::pair<int, modeMap>(index, result));
                        break;
                        }
                    case 2 :
                        {
                        std::lock_guard<std::mutex> guard(mapMutexThread2);
                        m_mapCollectThread2.insert(std::pair<int, modeMap>(index, result));
                        break;
                        }
                
                    case 3 : 
                        {
                        std::lock_guard<std::mutex> guard(mapMutexThread3);
                        m_mapCollectThread3.insert(std::pair<int, modeMap>(index, result));
                        break;
                        }
                    

                    case 4 :
                        {
                        std::lock_guard<std::mutex> guard(mapMutexThread4);
                        m_mapCollectThread4.insert(std::pair<int, modeMap>(index, result));
                        break;
                        }	
                    }	
                return;
            }
            
            
              /**
            
            Recursive function to sub matrix into 2x2 blocks which are then used to compute downsampled images.            
             
             @param 
                    rowStart, rowEnd - row positions of the elements in the baseImage 
                    colStart, colEnd - col psitiiop pf tje elemtnes in the baseImage.
                    depth - Depth of the cube in the baseImage 
                    threadNumber - thread reading a particular 2x2 block in the baseImage.

             @return 
                    void           
            */
            
            void divideCube(int rowStart, int rowEnd, int colStart, int colEnd, int depth, int threadNumber)
            {

                if((rowEnd-rowStart) == 2 && (colEnd-colStart) == 2)
                {	
                    modeMap result = findMode(rowStart, rowEnd, colStart, colEnd, depth, threadNumber);
                    int index = (rowEnd/2)*(m_dimB/2) + (colEnd/2);
                    addBlock(threadNumber, index, result);
                    return;
                } 

                else if((rowEnd - rowStart) !=2)
                {
                    divideCube( rowStart, rowStart + (rowEnd-rowStart)/2, colStart, colEnd, depth, threadNumber);


                    divideCube( rowStart + (rowEnd-rowStart)/2, rowEnd, colStart, colEnd, depth, threadNumber);
                
                }
                else if((colEnd-colStart) != 2) 
                {       
                          divideCube( rowStart, rowEnd, colStart, colStart + (colEnd - colStart)/2, depth, threadNumber);
                                
                          divideCube( rowStart, rowEnd, colStart+(colEnd - colStart)/2, colEnd, depth, threadNumber);
                
                }

                return;
            }
            
             /**
            
            Function to merge all m_mapCollectThread1(2,3,4) map objects and populate globalMap object.        
            Memory is freed from all the m_mapCollectThread1(2,3,4)  objects. 
             
             @param 
                    void  
             @return 
                    void           
            */
            
            void mergeAllMaps()
            {    


                m_globalMap.insert(m_mapCollectThread1.begin(), m_mapCollectThread1.end());
                m_mapCollectThread1.clear();
                
                m_globalMap.insert(m_mapCollectThread2.begin(), m_mapCollectThread2.end());
                m_mapCollectThread2.clear();
                
                m_globalMap.insert(m_mapCollectThread3.begin(), m_mapCollectThread3.end());
                m_mapCollectThread3.clear();
                
                m_globalMap.insert(m_mapCollectThread4.begin(), m_mapCollectThread4.end());
                m_mapCollectThread4.clear();
                
                //Free the memory allocated for the image since we now have all the data in globalMap
                m_baseImage.resize(extents[0][0]);
                
                return;
            }
            
            //Getter function to read globalMap size. 
            size_t getGlobalMapSize()
            {
                return m_globalMap.size();
            }
            
            //Getter function to read depth
            int getDepth()
            {
                return m_depth;
            }
            
            /**
            
            Function to print all downsampled versions of a 2D image. This function also calls the reduceGlobalMap() which groups 2x2 blocks from the base image.
            It runs through a loop and print values from the globalMap. It also keeps track of rows and colums of the downsampled images
             
             @param 
                    void  
             @return 
                    void           
            */
            void printDownsampled()
            {
                
                //Outerloop to print the number of downsampled images based on the value of minDim.
                for(int k = m_minDim; k > 1; k /= 2)
                {
                    int j = 1;
                    
                    //Inner loop to print elements 
                    for(auto i = m_globalMap.begin(); i != m_globalMap.end(); ++i, ++j)
                    {
                        
                        std::cout << (i->second).getMode() << " " ;
                        
                        if (j%(m_downCols) == 0 ) 
                            std::cout << std::endl;
                    }
                    
                    std::cout << std::endl;
                    
                    //Reduce the globalMap size after every iteration. This increases the grouping size in the base image and forms new Map objects. 
                    if(getGlobalMapSize() >1) reduceGlobalMap();
                    
                }
                
                 return;
            }
            
            
            /**
            
            Helper funtion for printDownsampled(). Everytime reduceGlobalMap() is called it groups 2x more elements in row and columns of the base image. 
            Since we are keep track of element counts in individual modeMaps, we need not read the base Image everytime to produce a downsampled image. 
             
             
             From the description of class twoDArray, this is a 1-downsampled 2x2 image where each element is a modeMap object having data about the different elements in the baseImage and their counts.
             
             01 02 | 03 04
             05 06 | 07 08
             -----------------
             09 10 | 11 12
             13 14 | 15 16
             
             After downsampling, the new set of modeMap objects look like this 
             
             01new | 02new
             -----------------
             03new | 04new
             
             where,

             01new = group ( 01, 02, 05, 06)
             02new = group ( 03, 04, 07, 08)
             03new = group ( 09, 10, 13, 14)
             04new = group ( 11, 12, 15, 16)
             
             And the process continues until one we reach the minDim.
             
             @param 
                    void  
             @return 
                    void           
            */
            
            void reduceGlobalMap()
            {
                int j =1;
                map < int, modeMap > new_globalMap;
                modeMap updatedCube;
                for(auto i = m_globalMap.begin(); i != m_globalMap.end(); ++j)
                {
                    
                    auto copy = i, k =i, l =i, m = i;
                    advance(k,1);
                    advance(l,m_downCols);
                    advance(m,m_downCols+1);
        
                    //Add for modeMaps to form a new one for the downsampled image
                    updatedCube = (i->second + k->second)+ (l->second + m->second) ;                  
                    
                    //Update the mode for the new modeMap
                    updatedCube.calculateMode();
                    
                    //Advance the current iterator by 4 modeMaps since we are grouping only 4 modeMaps at any step in the process.
                    advance(i,4);
        
                    // form a new globalMap 
                    new_globalMap.insert(std::pair<int, modeMap>(j, updatedCube));        
                    
                }
                
                //Clear memory once downsampling is done and update member variables accordingly. 
                m_globalMap.clear();
                m_globalMap = new_globalMap;
                m_downCols /= 2;
                m_downRows /= 2;
                m_minDim /= 2; 
               
               std::cout << "The number of cubes is : " << getGlobalMapSize()  << std::endl;               
               
               return;
                
            }
};

/**
            
            Recursive function to split a 2D array into 4 symmetric blocks and have one thread working on each of them. 
            Recursion bottoms out when max threads have been created or min block size of 2x2 is reached.
             
             @param 
                    rowStart, rowEnd - row positions of the elements in the baseImage 
                    colStart, colEnd - col psitiiop pf tje elemtnes in the baseImage.
                    depth - Depth of the cube in the baseImage 
                    threadNumber - thread reading a particular 2x2 block in the baseImage.  
                    userImage - a twoDArray object whose downsampled versions we are interested in.
                    totalThreads - totalThreads that have been created by this function
                    
             @return 
                    void           
            */
            
void  startThreading(twoDArray & userImage, int rowStart, int rowEnd, int colStart, int colEnd, int depth, int threadNumber = 1, int totalThreads = 1)
{

	if(depth > 1 && (rowEnd-rowStart) > 2 && (colEnd-colStart) > 2 && totalThreads < MAX_THREADS)
	{
 
		auto thread1 = std::thread(&startThreading, ref(userImage), rowStart, rowStart + (rowEnd - rowStart)/2, colStart, colStart + (colEnd - colStart)/2, depth-1, 1, totalThreads*4);

		auto thread2 = std::thread(&startThreading, ref(userImage), rowStart, rowStart + (rowEnd - rowStart)/2, colStart + (colEnd - colStart)/2, colEnd, depth-1, 2, totalThreads*4);

		auto thread3 = std::thread(&startThreading, ref(userImage), rowStart + (rowEnd - rowStart)/2, rowEnd, colStart, colStart + (colEnd - colStart)/2, depth-1, 3, totalThreads*4);

		auto thread4 = std::thread(&startThreading, ref(userImage), rowStart + (rowEnd - rowStart)/2, rowEnd, colStart + (colEnd - colStart)/2, colEnd, depth-1, 4, totalThreads*4);

		thread1.join();
        thread2.join();
		thread3.join();
		thread4.join();
		
	} 
	
	else
	{
		userImage.divideCube(rowStart, rowEnd, colStart, colEnd, depth, threadNumber);
	}
	return;	
}



int main()

{

// Timers to measure performance;  
high_resolution_clock::time_point t1 = high_resolution_clock::now();   

//Set the seed for random number generation 
srand(3);


int dimA, dimB ;


std::cout << "Enter the number of rows (in powers of 2) : " << std::endl;
std::cin >> dimA;

std::cout << "Enter the number of cols (in powers of 2) : "  << std::endl;
std::cin >> dimB;

//Create an object
twoDArray ImageData(dimA,dimB);

//Start populating the m_mapCollectThread1(2,3,4)  objects. 
startThreading(ImageData, 0, dimA, 0, dimB, ImageData.getDepth());

//Merge all m_mapCollectThread1(2,3,4)  objects together.
ImageData.mergeAllMaps();


//std::cout << "The number of smallest cubes is : " << ImageData.getGlobalMapSize()  << std::endl;
ImageData.printDownsampled();

high_resolution_clock::time_point t2 = high_resolution_clock::now();
auto duration = duration_cast<milliseconds>( t2 - t1 ).count();
std::cout << "The duration is : " << duration << " milli seconds " << std::endl;

return 0;

}

