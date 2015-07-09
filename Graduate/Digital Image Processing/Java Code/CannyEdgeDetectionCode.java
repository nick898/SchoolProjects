/**
 Name: Nick Kallfa
 Course: Digital Image Processing Fall 2014

 Description: This plug in filter implements the Canny Edge Detector using ImageJ. The code starts with two sections that contain preliminary steps before beginning to
              implement the Canny Edge Detector. These preliminary steps create a number of blank ImageProcessor objects to be used throughout the algorithm and the other preliminary
              step is to obtain user input to determine the size of the Gaussian kernel, the standard deviation, as well as the high and low threshold to be used in the last step
              of the Canny Edge Detection algorithm.

              The Canny Edge Detector is broken up into 4 steps:

              Step 1: Noise Reduction Via Gaussian
              Step 2: Compute Gradient Magnitude Image and Gradient Direction Image using Sobel Filters in the 0, 45, 90, and 135 degree directions
              Step 3: Non-Maximal Suppression
              Step 4: Thresholding with Hysteresis
 **/

import ij.IJ;
import ij.ImagePlus;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;
import ij.process.ByteProcessor;
import ij.process.FloatProcessor;
import ij.plugin.filter.Convolver;
import ij.gui.GenericDialog;
import ij.gui.MessageDialog;
import java.lang.Math;




public class HW_7 implements PlugInFilter {

    ImagePlus im; //New instance variable of this plugin object

	public int setup(String arg, ImagePlus im) {
	    if (im == null){
            IJ.noImage(); //currently no image is open
            return DONE;
	    }
	    this.im = im; //keep a reference to the image im
		return DOES_8G; // this plugin accepts 8-bit grayscale images
	}

	public void run(ImageProcessor ip) {
		int w = ip.getWidth(); //Get width of image
		int h = ip.getHeight(); //Get height of image


/**------------------------------------------------------------------  BEGIN PRELIMINARY STEPS ------------------------------------------------------------------**/

/**-----------------------------------  BEGIN: CREATE IMAGES THAT WILL BE USED IN COMPUTATIONS -----------------------------------**/
        // Create the smoothed image to be used for the different edge images
        ImageProcessor Smoothed_Ip_xf = new FloatProcessor(w, h);
        ImageProcessor Smoothed_Ip_yf = new FloatProcessor(w, h);
        ImageProcessor Smoothed_Ip_45f = new FloatProcessor(w, h);
        ImageProcessor Smoothed_Ip_135f = new FloatProcessor(w, h);


        // Create the edge image detecting vertical edges
        ImageProcessor G_xf_Ip = new FloatProcessor(w, h);

        // Create the edge image detecting horizontal edges
        ImageProcessor G_yf_Ip = new FloatProcessor(w, h);


        // Create the edge image detecting 135 degree edges
        ImageProcessor G_45f_Ip = new FloatProcessor(w, h);


        // Create the edge image detecting 45 degree edges
        ImageProcessor G_135f_Ip = new FloatProcessor(w, h);


        // Create the gradient magnitude image
        ImageProcessor GMag_Ip = new ByteProcessor(w, h); //Byte version
        ImageProcessor GMagf_Ip = new FloatProcessor(w, h); //Floating point version

        // Create the gradient direction image
        ImageProcessor GDir_Ip = new ByteProcessor(w, h); //Byte version
        ImageProcessor GDirf_Ip = new FloatProcessor(w, h); //Floating point version

        // Create the edge image (output from non-maximal suppression)
        ImageProcessor Edge_Ip = new ByteProcessor(w, h); //Byte Version
        Edge_Ip.setValue(0); // 0 = Black
        Edge_Ip.fill(); // Fill Edge_Ip all black

        ImageProcessor Copy_Edge_Ip = new ByteProcessor(w, h); //Byte Version
        Copy_Edge_Ip.setValue(0); // 0 = Black
        Copy_Edge_Ip.fill(); // Fill Edge_Ip all black

        ImageProcessor Edgef_Ip = new FloatProcessor(w, h); //Floating Point version
        Edgef_Ip.setValue(0); // 0 = Black
        Edgef_Ip.fill(); // Fill Edge_Ip all black

        ImageProcessor Copy_Edgef_Ip = new FloatProcessor(w, h); //Floating Point version
        Copy_Edgef_Ip.setValue(0); // 0 = Black
        Copy_Edgef_Ip.fill(); // Fill Edge_Ip all black

        // Create the Threshold with Hysteresis image
        ImageProcessor Threshold_Ip = new ByteProcessor(w, h);
        Threshold_Ip.setValue(0); // 0 = Black
        Threshold_Ip.fill(); // Fill Edge_Ip all black
/**-----------------------------------  END: CREATE IMAGES THAT WILL BE USED IN COMPUTATIONS -----------------------------------**/

/**-----------------------------------------------------  BEGIN: USER INPUT -----------------------------------------------------**/
        double STDDev = 0, TLow = 0, THigh = 0; //Declare and initialize variables for the standard deviation, low threshold, and high threshold
        int Size = 0; //Initialize size of Gaussian Filter
        boolean EdgeStrengthImage; // specifies whether edge strength image should be shown
        GenericDialog gd = new GenericDialog("User Inputs");
        gd.addNumericField("Size of Gaussian Filter (Odd Integer)", Size, 0); //Field for Size of Gaussian Filter
        gd.addNumericField("Standard Deviation", STDDev, 0); //Field for Standard Deviation
        gd.addNumericField("Low Threshold (1 - 255)", TLow, 0); //Field for Low Threshold
        gd.addNumericField("High Threshold (1 - 255)", THigh, 0); //Field for High Threshold
		gd.showDialog();
		if (gd.wasCanceled())
		{
			return;
		}
		else
		{
            Size = (int) gd.getNextNumber(); //Set Size variable from user input. This allows the user to set the size of the Gaussian Filter
			STDDev = gd.getNextNumber(); //Set STDDev variable from user input
            TLow = gd.getNextNumber(); //Set TLow variable from user input
			THigh = gd.getNextNumber(); //Set THigh variable from user input
		}

/**-----------------------------------------------------  END: USER INPUT -----------------------------------------------------**/

/**------------------------------------------------------------------  END PRELIMINARY STEPS ------------------------------------------------------------------**/

/**------------------------------------------------------------------  BEGIN CANNY EDGE DETECTION ------------------------------------------------------------------**/

/**-----------------------------------------------------   BEGIN STEP 1: NOISE REDUCTION VIA GAUSSIAN ----------------------------------------------------- **/

        ImageProcessor ipf = ip.convertToFloat(); //Convert original image to floating point
        int SizeSquared = (int) Math.pow(Size,2); //Square the size of Gaussian Kernel
        int HalfSize = (Size - 1)/2; //Cut the size of the Gaussian Kernel almost in half
        float pix; //Temporary storage to be used throughout code to store pixel values
        float[] GaussianFilter = new float[SizeSquared]; //Initialize Gaussian Filter to be used in the convolution
        double[] GaussianFilterD = new double[SizeSquared]; //Initialize Gaussian Filter to be used...this filter will be composed of numbers of type double
        double Constant = 1/(2*Math.PI*Math.pow(STDDev,2)); //Compute 1/(2*pi*sigma^2)
        double ExponentDenom = 2*Math.pow(STDDev,2); //Compute 2*sigma^2
        double Value; //Temporary storage to store the computations that will form the Gaussian Filter

        //FOR LOOP TO FORM GaussianFilterD
        for(int i = 0; i < Size; i++)
        {
            for(int j = 0; j < Size; j++)
            {
                Value = Math.exp(-1*(Math.pow(j-HalfSize,2) + Math.pow(i-HalfSize,2))/(ExponentDenom)); //Set Value = e^(-(i^2 + j^2)/(2*sigma^2))
                GaussianFilterD[Size*i+j] = Constant*Value; //Place Value in GaussianFilterD
            }
        }

        //FOR LOOP TO FORM GaussianFilter using GaussianFilterD
        for(int i = 0; i < Math.pow(Size,2); i++)
        {
            GaussianFilter[i] = (float) GaussianFilterD[i]; //Convert double values to float one by one
        }

        //CONVOLVE IMAGE WITH GAUSSIAN
        Convolver cv = new Convolver(); //Create the convolver
        cv.setNormalize(true); //Normalize the filter
        cv.convolve(ipf,GaussianFilter,Size,Size); //Apply the GaussianFilter using convolution on the image ipf

        //For loop to create 4 smoothed images to be used in detecting the 4 edges: 0 deg, 45 deg, 90 deg, 135 deg
        for(int i = 0; i < w; i++)
        {
            for(int j = 0; j < h; j++)
            {
                pix = ipf.getf(i,j);
                Smoothed_Ip_xf.setf(i,j,pix); //Smoothed image to convolve with Sobel filter in x-direction (detects vertical edges)
                Smoothed_Ip_yf.setf(i,j,pix); //Smoothed image to convolve with Sobel filter in y-direction (detects horizontal edges)
                Smoothed_Ip_45f.setf(i,j,pix); //Smoothed image to convolve with Sobel filter in 45 degree direction (detects 135 degree edges)
                Smoothed_Ip_135f.setf(i,j,pix); //Smoothed image to convolve with Sobel filter in 135 degree direction (detects 45 degree edges)

            }
        }

/**-----------------------------------------------------   END STEP 1: NOISE REDUCTION VIA GAUSSIAN ----------------------------------------------------- **/

/**-----------------------------------------------------   BEGIN STEP 2: COMPUTE GRADIENT MAGNITUDE AND DIRECTION IMAGES ----------------------------------------------------- **/

        //Sobel Filters to detect edges
        float[] G_x = {-1,0,1,-2,0,2,-1,0,1}; //Sobel operator in x direction
        float[] G_y = {-1, -2, -1,0,0,0,1,2,1}; //Sobel operator in y direction
        float[] G_135 = {2,1,0,1,0,-1,0,-1,-2}; //Sobel operator in 135 degree direction
        float[] G_45 = {0, 1, 2,-1,0,1,-2,-1,0}; //Sobel operator in 45 degree direction

        //CONVOLVE IMAGE USING SOBEL FILTERS
        cv.convolve(Smoothed_Ip_xf,G_x,3,3); //Apply the Sobel filter in x-direction using convolution on the smoothed image
        cv.convolve(Smoothed_Ip_yf,G_y,3,3); //Apply the Sobel filter in y-direction using convolution on the smoothed image
        cv.convolve(Smoothed_Ip_45f,G_45,3,3); //Apply the Sobel filter in 45 degree direction using convolution on the smoothed image
        cv.convolve(Smoothed_Ip_135f,G_135,3,3); //Apply the Sobel filter in 135 degree direction using convolution on the smoothed image


        //For loop to define the 4 floating point images G_xf_Ip, G_yf_Ip, G_45f_Ip, G_135f_Ip
        for(int i = 0; i < w; i++)
        {
            for(int j = 0; j < h; j++)
            {
                pix = Smoothed_Ip_xf.getf(i,j);
                G_xf_Ip.setf(i,j,pix);
                pix = Smoothed_Ip_yf.getf(i,j);
                G_yf_Ip.setf(i,j,pix);
                pix = Smoothed_Ip_45f.getf(i,j);
                G_45f_Ip.setf(i,j,pix);
                pix = Smoothed_Ip_135f.getf(i,j);
                G_135f_Ip.setf(i,j,pix);
            }
        }

        //COMPUTE THE GRADIENT MAGNITUDE IMAGE GMagf_Ip
        float pix1, pix2; //Variables to store pixel values of G_xf_Ip and G_yf_Ip
        for(int i = 0; i < w; i++)
        {
            for(int j = 0; j < h; j++)
            {
                pix1 = G_xf_Ip.getf(i,j); //Get pixel value
                pix2 = G_yf_Ip.getf(i,j); //Get Pixel Value
                pix  = (float) Math.sqrt(pix1*pix1 + pix2*pix2); //Take the square root of the sum of the squares of both pixel values
                GMagf_Ip.setf(i,j,255*pix); //Place in (i,j) position of GMagf_Ip and scale by 255 to view
            }

        }


        //COMPUTE THE GRADIENT DIRECTION IMAGE
        double v1, v2, v3, v4; //Temporary storage for the abs. value of pixel vlaues
        int pixi; //Temporary storage for pixel values
        double[][] GDir = new double[h][w]; //Define 2D Array to store gradient direction values. Will use this during Step 3 Non-Maximal Suppression instead of the gradient direction image

        //For loop to compute gradient direction image
        for(int i = 0; i < w; i++)
        {
            for(int j = 0; j < h; j++)
            {
                //Obtain the absolute value of each pixel in the following four floating point images
                v1 = Math.abs(G_xf_Ip.getf(i,j));
                v2 = Math.abs(G_yf_Ip.getf(i,j));
                v3 = Math.abs(G_45f_Ip.getf(i,j));
                v4 = Math.abs(G_135f_Ip.getf(i,j));

                //If statements to find the maximal response (strongest edge direction)

                if (v1 > v2 && v1 > v3 && v1 > v4) //Vertical edge orientation
                {
                    GDir_Ip.putPixel(i,j,0);
                    GDir[j][i] = 0; //0 stands for vertical edge orientation
                }
                else if (v3 > v1 && v3 > v2 && v3 > v4) //45 degree edge orientation
                {
                    GDir_Ip.putPixel(i,j,1);
                    GDir[j][i] = 1; //1 stands for 45 degree edge orientation
                }
                else if (v2 > v1 && v2 > v3 && v2 > v4) //Horizontal edge orientation
                {
                    GDir_Ip.putPixel(i,j,2);
                    GDir[j][i] = 2; //2 stands for horizontal edge orientation
                }
                else if (v4 > v1 && v4 > v2 && v4 > v3) //135 degree edge orientation
                {
                    GDir_Ip.putPixel(i,j,3);
                    GDir[j][i] = 3; //3 stands for 135 degree edge orientation
                }

                pixi = (int) (255.0/3.0)*GDir_Ip.getPixel(i,j); //Scale the gradient direction image so we can actually see it
                GDir_Ip.putPixel(i,j,pixi); //Place pixi in the (i,j) location of GDir_Ip...we can now view the gradient direction image
            }
        }

                //OUTPUT GRADIENT DIRECTION IMAGE
                String DirTitle = "Gradient Direction";
				ImagePlus GDir_Im = new ImagePlus(DirTitle, GDir_Ip);
                GDir_Im.show();


                //OUTPUT GRADIENT MAGNITUDE IMAGE (need to convert to ByteProcessor to view first)
                GMagf_Ip.resetMinAndMax();
                GMag_Ip.insert(GMagf_Ip.convertToByte(true), 0, 0);

                String MagTitle = "Gradient Magnitude";
				ImagePlus GMag_Im = new ImagePlus(MagTitle, GMag_Ip);
                GMag_Im.show();

/**-----------------------------------------------------   BEGIN STEP 2: COMPUTE GRADIENT MAGNITUDE AND DIRECTION IMAGES ----------------------------------------------------- **/

/**-----------------------------------------------------  BEGIN STEP 3: NON-MAXIMAL SUPPRESSION -----------------------------------------------------**/
        // Ignoring boundaries for convenience
        float Magnitude; //Storage for magnitude
        double Direction; //Storage for direction

        //FOR LOOP TO COMPUTE NON-MAXIMAL SUPPRESSION IMAGE
        for(int i = 1; i < w - 1; i++)
        {
            for(int j = 1; j < h - 1; j++)
            {
                Magnitude = GMagf_Ip.getf(i,j); //Get the magnitude in the (i,j) position of the gradient magnitude image

                if (Magnitude != 0 ) //If the magnitude is non-zero then we find the direction of the edge
                {
                    Direction = GDir[j][i]; //Obtain the direction
                    float n1GradMag = 0, n2GradMag = 0; //Initialize storage for the magnitude of the 2 neighboring pixels

                    if(Direction == 0) //If Direction is 0 get the magnitude in the columns to the left and right
                    {
                        n1GradMag = GMagf_Ip.getf(i-1,j);
                        n2GradMag = GMagf_Ip.getf(i+1,j);
                    }
                    else if (Direction == 1) //If Direction is 45 degree get the magnitude in adjacent pixels
                    {
                        n1GradMag = GMagf_Ip.getf(i+1,j-1);
                        n2GradMag = GMagf_Ip.getf(i-1,j+1);
                    }
                    else if (Direction == 2) //If Direction is 2 get the magnitude in the rows above and below
                    {
                        n1GradMag = GMagf_Ip.getf(i,j-1);
                        n2GradMag = GMagf_Ip.getf(i,j+1);
                    }
                    else if (Direction == 3) //If Direction is 135 degrees get the magnitude in adjacent pixels
                    {
                        n1GradMag = GMagf_Ip.getf(i-1,j-1);
                        n2GradMag = GMagf_Ip.getf(i+1,j+1);
                    }

                    if (Magnitude > n1GradMag && Magnitude > n2GradMag) //Check to see if the magnitude of pixel under inspection is the largest relative to its adjacent pixels
                    {
                        pix = GMagf_Ip.getf(i,j); //Store magnitude in (i,j) position in pix variable
                        Edgef_Ip.setf(i,j,pix); //Place this pixel value in the (i,j) position of Edge Image (Edge Image is output from non-maximal suppression)
                        Copy_Edgef_Ip.setf(i,j,pix); //Place this pixel value in the (i,j) position of Copy Edge Image (this will be used for thresholding with hysteresis)

                    }
                    else{} //Else do nothing since edge image was already filled with 0's at beginning of code
                }
                else{} //Do nothing if magnitude is 0
            }
        }


                Copy_Edgef_Ip.resetMinAndMax();
                Copy_Edge_Ip.insert(Copy_Edgef_Ip.convertToByte(true), 0, 0);

                //Display edge magnitude image but we need to convert to ByteProcessor first
                Edgef_Ip.resetMinAndMax();
                Edge_Ip.insert(Edgef_Ip.convertToByte(true), 0, 0);

                String cTitle = "Non-Maximal Suppression";
				ImagePlus Edge_Im = new ImagePlus(cTitle, Edge_Ip);
                Edge_Im.show();
/**-----------------------------------------------------  END STEP 3: NON-MAXIMAL SUPPRESSION -----------------------------------------------------**/

/**-----------------------------------------------------  BEGIN STEP 4: THRESHOLDING WITH HYSTERESIS -----------------------------------------------------**/
        int count = 0; //Initialize counter to be used in while loop
        int IterationCount = 500; //Number of iterations to perform

        //Scan through all pixels in the edge image and mark all edges with magnitude above the high threshold as a true edge otherwise if the magnitude is below the low threshold
        //then we delete that pixel (set it to zero)
        for(int i = 0; i < w; i++)
        {
            for(int j = 0; j < h; j++)
            {
                Magnitude = Copy_Edge_Ip.getPixel(i,j); //Obtain magnitude in edge image

                if (Magnitude > THigh) //If Magnitude is larger than the high threshold then make pixel white
                {
                    Threshold_Ip.putPixel(i,j,255); //True edge (Updates Threshold image to be output)
                    Copy_Edge_Ip.putPixel(i,j,255); //True edge (Updates edge image which will be used in the iterative thresholding)
                }
                else if(Magnitude < TLow) //Else if magnitude is below the low threshold then make pixel black
                {
                    Threshold_Ip.putPixel(i,j,0); //Not an edge (Updates Threshold image to be output)
                    Copy_Edge_Ip.putPixel(i,j,0); //Not an edge (Updates edge image which will be used in the iterative thresholding)
                }
            }
        }

    while(count < IterationCount) //Iterate again and again
    {
        //Ignore boundary pixels for convenience
        for(int i = 1; i < w - 1; i++)
        {
            for(int j = 1; j < h - 1; j++)
            {
                Magnitude = Copy_Edge_Ip.getPixel(i,j); //Obtain magnitude in edge image
                if (Magnitude == 255) //If we reach a true edge then we look at its 8 neighbors
                {
                    //Obtain the magnitude in the 8 neighbors
                    int n1,n2,n3,n4,n5,n6,n7,n8;
                    n1 = Copy_Edge_Ip.getPixel(i-1,j);
                    n2 = Copy_Edge_Ip.getPixel(i-1,j-1);
                    n3 = Copy_Edge_Ip.getPixel(i,j-1);
                    n4 = Copy_Edge_Ip.getPixel(i+1,j-1);
                    n5 = Copy_Edge_Ip.getPixel(i+1,j);
                    n6 = Copy_Edge_Ip.getPixel(i+1,j+1);
                    n7 = Copy_Edge_Ip.getPixel(i,j+1);
                    n8 = Copy_Edge_Ip.getPixel(i-1,j+1);
                    if(n1 >= TLow) //If n1 is greater than or equal to the low threshold then we mark it as an edge
                    {
                        Copy_Edge_Ip.putPixel(i-1,j,255); //Update edge image
                        Threshold_Ip.putPixel(i-1,j,255); //Update threshold image (This is the output image from thresholding with hysteresis)
                    }
                    if(n2 >= TLow) //If n2 is greater than or equal to the low threshold then we mark it as an edge
                    {
                        Copy_Edge_Ip.putPixel(i-1,j-1,255); //Update edge image
                        Threshold_Ip.putPixel(i-1,j-1,255); //Update threshold image (This is the output image from thresholding with hysteresis)
                    }
                    if(n3 >= TLow) //If n3 is greater than or equal to the low threshold then we mark it as an edge
                    {
                        Copy_Edge_Ip.putPixel(i,j-1,255); //Update edge image
                        Threshold_Ip.putPixel(i,j-1,255); //Update threshold image (This is the output image from thresholding with hysteresis)
                    }
                    if(n4 >= TLow) //If n4 is greater than or equal to the low threshold then we mark it as an edge
                    {
                        Copy_Edge_Ip.putPixel(i+1,j-1,255); //Update edge image
                        Threshold_Ip.putPixel(i+1,j-1,255); //Update threshold image (This is the output image from thresholding with hysteresis)
                    }
                    if(n5 >= TLow) //If n5 is greater than or equal to the low threshold then we mark it as an edge
                    {
                        Copy_Edge_Ip.putPixel(i+1,j,255); //Update edge image
                        Threshold_Ip.putPixel(i+1,j,255); //Update threshold image (This is the output image from thresholding with hysteresis)
                    }
                    if(n6 >= TLow) //If n6 is greater than or equal to the low threshold then we mark it as an edge
                    {
                        Copy_Edge_Ip.putPixel(i+1,j+1,255); //Update edge image
                        Threshold_Ip.putPixel(i+1,j+1,255); //Update threshold image (This is the output image from thresholding with hysteresis)
                    }
                    if(n7 >= TLow) //If n7 is greater than or equal to the low threshold then we mark it as an edge
                    {
                        Copy_Edge_Ip.putPixel(i,j+1,255); //Update edge image
                        Threshold_Ip.putPixel(i,j+1,255); //Update threshold image (This is the output image from thresholding with hysteresis)
                    }
                    if(n8 >= TLow) //If n8 is greater than or equal to the low threshold then we mark it as an edge
                    {
                        Copy_Edge_Ip.putPixel(i-1,j+1,255); //Update edge image
                        Threshold_Ip.putPixel(i-1,j+1,255); //Update threshold image (This is the output image from thresholding with hysteresis)
                    }
                }
            }
        }
        count++; //Update counter and continue iterations
    }

    //Display the threshold with hysteresis image.
    String ThreshTitle = "Threshold with Hysteresis";
    ImagePlus Threshold_Im = new ImagePlus(ThreshTitle, Threshold_Ip);
    Threshold_Im.show();
/**-----------------------------------------------------  END STEP 4: THRESHOLDING WITH HYSTERESIS -----------------------------------------------------**/

/**------------------------------------------------------------------  END CANNY EDGE DETECTION ------------------------------------------------------------------**/

}
}
