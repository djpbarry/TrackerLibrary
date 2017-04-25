package ProbabilisticTracking;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.GenericDialog;
import ij.gui.StackWindow;
import ij.plugin.ZProjector;
import ij.plugin.filter.PlugInFilter;
import ij.process.ByteProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import ij.text.TextWindow;
import ij.gui.ImageCanvas;
import ij.io.FileInfo;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Vector;
import java.util.Random;
import java.util.regex.Pattern;
import java.awt.Button;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.GridLayout;
import java.awt.Panel;
import java.awt.Point;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;

/**
 * <h2>PFTracking3D</h2> <h3>An ImageJ Plugin base class for tracking in 4D
 * fluorescence microscopy</h3>
 *
 * <p>The provided algorithm serves to model complex models with (hopefully) a
 * lot of a-priori knowledge of the dynamics. If you're able to crop a object to
 * track, you might use this algorithm with a random walk model with large
 * variances in the dynamics distribution. The algorithm is able to handle very
 * low signal to noise ratios (a feature point can be tracked in a 2D image till
 * a SNR of 2.5; for 3D images the SNR might be even smaller). If you have many
 * objects close to each other doing large and unknown movements, this is not
 * the right choice.
 *
 * <p>Correctly overwriting the abstract methods should track any object.
 * Basically, the dimensions of your state vector have to be set, and from there
 * on, a expected image from such a state have to be created in the
 * <code>generateIdealImage</code> method. Further you should define the
 * proposal distribution by overwriting the method
 * <code>drawFromProposalDistribution</code>. As an example of implementation,
 * have a look on the
 * <code>LinearMovementFPTracker_3D</code> class or FPTracker_3D class.
 *
 * <p>The data structures provide support for multiple objects. Note that multi
 * target tracking (MTT) only does not work properly with this algorithm if the
 * objects are too near from each other.
 *
 * <p>For further informations please consider the appropriate tutorial.
 *
 * <p> The tracks can be processed on 8-bit, 16-bit and 32-bit but no color
 * images.
 *
 * <p><b>Disclaimer</b> <br>IN NO EVENT SHALL THE ETH BE LIABLE TO ANY PARTY FOR
 * DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING
 * LOST PROFITS, ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION,
 * EVEN IF THE ETH HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. THE ETH
 * SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
 * THE SOFTWARE PROVIDED HEREUNDER IS ON AN "AS IS" BASIS, AND THE ETH HAS NO
 * OBLIGATIONS TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR
 * MODIFICATIONS.<p>
 *
 * @version 1.0. January, 2008 (requires: Java 5 or higher)
 * @author Janick Cardinale - Master student at the
 * <a href="http://www.cbl.ethz.ch/">Computational Biophysics Lab</a>, ETH Zurich
 */
public abstract class PFTracking3D implements PlugInFilter {

    protected static enum STATE_OF_FILTER {

        WAITING, INIT, READY_TO_RUN, RUNNING, VISUALIZING, CORRECTING
    };
    protected static final String RESULT_FILE_SUFFIX = "_mtTracker_results.txt";
    protected static final String INIT_FILE_SUFFIX = "_mtTracker_initValues.txt";
    protected STATE_OF_FILTER mStateOfFilter = STATE_OF_FILTER.WAITING;
    /*
     * Parameters
     */
    protected int mNbThreads = Runtime.getRuntime().availableProcessors();
//	protected int mNObjects = 1;
    protected int mNbParticles = 1500;
    protected int mRepSteps = 5;
    protected int mInitRWIterations = 1;
    protected int mResamplingThreshold = mNbParticles / 2;
    protected float mBackground = 1f;
    protected float[] mSigmaOfRandomWalk;
    protected String[] mDimensionsDescription;
    protected float mSigmaPSFxy = 192f;
    protected float mSigmaPSFz = 268f;//in nm
    protected long mSeed = 88888888;
    protected int mWavelengthInNm = 450;
    protected float mNA = 1.2f;
    protected float mn = 1.3f;
    protected int mTrackTillFrameNb = 0;

    /*
     * Options
     */
    protected boolean mDoResampling = true;
    protected boolean mDoPrintStates = true;
    protected boolean mDoPrecisionOptimization = true;

    /*
     * Monitoring variables
     */
    protected boolean mDoMonitorIdealImage = false;
    protected boolean mDoMonitorParticles = false;
    protected ImageStack mIdealImageMonitorStack;
    /**
     * Frames, Objects, Particles, Particle(x,y,Intensity,weight)
     */
    protected Vector<Vector<Vector<float[]>>> mParticleMonitor = new Vector<Vector<Vector<float[]>>>(); //init grosse?	
    //
    // frequently used members(for simplicity)
    //
    /**
     * Starts with 1..NFrames
     */
    protected int mFrameOfInitialization = 1;
    protected int mHeight, mWidth, mNSlices, mNFrames;
    protected float mPxWidthInNm, mPxDepthInNm;
    protected Random mRandomGenerator = new Random(mSeed);
    protected ImagePlus mOriginalImagePlus;
    protected ImagePlus mZProjectedImagePlus;
    protected Vector<float[]> mStateVectors = new Vector<float[]>();
    /**
     * Stores the state vectors for each frame
     */
    protected Vector<Vector<float[]>> mStateVectorsMemory = new Vector<Vector<float[]>>(); //grosse->evt in setup() initialisieren
    protected float[] mMaxLogLikelihood;
    protected Vector<Vector<float[]>> mParticles = new Vector<Vector<float[]>>();

    /**
     * The plugins setup method, invoked by imageJ
     */
    public int setup(String aArgs, ImagePlus aImp) {
        if (IJ.versionLessThan("1.38u")) {
            return DONE;
        }
        if (aImp == null) {
            IJ.showMessage("Please open an image to track first.");
            return DONE;
        }

        if (aImp.getNFrames() < 2) {
            //nothing to track
//			IJ.showMessage("The image only contains one frame, check your image properties");
            IJ.run("Properties...");
//			return DONE;
        }

        while (true) {
            String vUnit = aImp.getCalibration().getUnit();
            if (vUnit.equals("nm")) {
                mPxWidthInNm = (float) aImp.getCalibration().pixelWidth;
                mPxDepthInNm = (float) aImp.getCalibration().pixelDepth;
                break;
            } else if (vUnit.equals(IJ.micronSymbol + "m")) {
                mPxWidthInNm = (float) aImp.getCalibration().pixelWidth * 1000;
                mPxDepthInNm = (float) aImp.getCalibration().pixelDepth * 1000;
                break;
            } else if (vUnit.equals("mm")) {
                mPxWidthInNm = (float) aImp.getCalibration().pixelWidth * 1000000;
                mPxDepthInNm = (float) aImp.getCalibration().pixelDepth * 1000000;
                break;
            }
            IJ.showMessage("Please enter the pixel sizes in nm, " + IJ.micronSymbol + "m or mm");
            IJ.run("Properties...");
        }

        //
        // Init members
        //
        mOriginalImagePlus = aImp;
        mHeight = mOriginalImagePlus.getHeight();
        mWidth = mOriginalImagePlus.getWidth();
        mNFrames = mOriginalImagePlus.getNFrames();
        mNSlices = mOriginalImagePlus.getNSlices();
        mTrackTillFrameNb = mNFrames;
        mStateVectorsMemory.setSize(mOriginalImagePlus.getNFrames());
        mMaxLogLikelihood = new float[mOriginalImagePlus.getNFrames()];
        mFrameOfInitialization = sliceToFrame(mOriginalImagePlus.getCurrentSlice());
        mRandomGenerator = new Random(mSeed);

        mSigmaPSFxy = (0.21f * mWavelengthInNm / mNA);
        mSigmaPSFz = (0.66f * mWavelengthInNm * mn / (mNA * mNA));

        mDimensionsDescription = getMDimensionsDescription();
        mSigmaOfRandomWalk = getMSigmaOfRandomWalk();

        if (!getUserDefinedParams()) {
            return DONE;
        }
        doZProjection();
        initMonitoring();
        initVisualization();

        initPlugin();

        if (mStateOfFilter != STATE_OF_FILTER.READY_TO_RUN) {
            return DONE;
        }
        return DOES_8G + DOES_16 + DOES_32 + STACK_REQUIRED + PARALLELIZE_STACKS;
    }

    public void run(ImageProcessor ip) {
        if (mStateOfFilter != STATE_OF_FILTER.RUNNING
                && //mStateOfFilter != STATE_OF_FILTER.VISUALIZING && 
                mStateOfFilter != STATE_OF_FILTER.READY_TO_RUN) {
            IJ.showMessage("No valid initialization. No calculation started");
            return;
        }

        if (!showParameterDialog()) {
            return;
        }

        initParticleFilter(getAFrameCopy(mOriginalImagePlus, mFrameOfInitialization), mInitRWIterations, mFrameOfInitialization);

        initVisualization();
        mStateVectorsMemory.setElementAt(copyStateVector(mStateVectors), mFrameOfInitialization - 1);

        runParticleFilter(mOriginalImagePlus);

        //
        // save or write data to screen
        //
        if (mDoPrintStates) {
            if (mOriginalImagePlus.getOriginalFileInfo() == null) {
                printStatesToWindow("Tracking states", mStateVectorsMemory);
            } else {
                writeResultFile(getResultFile());
            }
        }

        visualizeMonitors();
    }

    /**
     * The method first checks if there are already results found in the same
     * directory as the movie is saved. If yes, the visualizing mode is used. If
     * not it is checked if there is a initialization file in this directory. If
     * yes, the tracker automatically begins; this may be used to first
     * initialize several files and then let them track in a macro. If no
     * initialization was found, the method tries to automatically initialize
     * using a method one may override
     * <code>autoInitFilter</code>.
     */
    private void initPlugin() {
        //
        // Check for a result file
        //
        FileInfo vFI = mOriginalImagePlus.getOriginalFileInfo();
        if (vFI != null && !(vFI.directory == "")) {//|| vFI.fileName == "")) {
            //
            // Check if there are results to visualize
            //
            File vResFile = getResultFile();

            if (vResFile.exists()) {
                //read out all the stuff; parameter too? no
                if (readResultFile(vResFile)) {
                    mStateOfFilter = STATE_OF_FILTER.VISUALIZING;
                    return;
                }
            }

            //
            // Check if there are init values
            ///
            File vInitFile = getInitFile();

            if (vInitFile.exists()) {
                //read out all the stuff; parameter too? no
                if (readInitFile(vInitFile)) {
                    mStateOfFilter = STATE_OF_FILTER.READY_TO_RUN;
                    return;
                }
            }
        }
        //
        // try to init automatically
        //
        ImageStack vInitStack = getAFrameCopy(mOriginalImagePlus, mFrameOfInitialization);
        if (autoInitFilter(vInitStack)) {
            /*
             * Save the state vector after initialisation(the first frame is
             * finished)
             */
            mStateVectorsMemory.set(mFrameOfInitialization - 1, copyStateVector(mStateVectors));
            //we propose the initialisation but do not run
            mStateOfFilter = STATE_OF_FILTER.VISUALIZING;

            mZProjectedImagePlus.setSlice(sliceToFrame(mOriginalImagePlus.getCurrentSlice()));
            mZProjectedImagePlus.repaintWindow();
            return;
        }
        mStateOfFilter = STATE_OF_FILTER.WAITING;
    }

    private void initVisualization() {
        // generate the previewCanvas - while generating it the drawing will be done 
        DrawCanvas vDrawCanvas = new DrawCanvas(mZProjectedImagePlus);

        // display the image and canvas in a stackWindow  
        new TrajectoryStackWindow(mZProjectedImagePlus, vDrawCanvas);
    }

    private void initParticleFilter(ImageStack aInitStack, int aInitParticleFilterIterations, int aFrameOfInit) {
        //
        // - set up state vector
        // - create particles
        // - filter the initialized values
        //
        createParticles(mStateVectors, mParticles);
        filterTheInitialization(aInitStack, aInitParticleFilterIterations, aFrameOfInit);
    }

    /**
     * To override.
     *
     * @return true if successful.
     */
    protected boolean autoInitFilter(ImageStack aImageStack) {
        return false;
    }

    private void doZProjection() {
        ImageStack vZProjectedStack = new ImageStack(mWidth, mHeight);
        ZProjector vZProjector = new ZProjector(mOriginalImagePlus);
        vZProjector.setMethod(ZProjector.MAX_METHOD);
        for (int vC = 0; vC < mOriginalImagePlus.getNFrames(); vC++) {
            vZProjector.setStartSlice(vC * mNSlices + 1);
            vZProjector.setStopSlice((vC + 1) * mNSlices);
            vZProjector.doProjection();
            vZProjectedStack.addSlice("", vZProjector.getProjection().getProcessor());
        }
        mZProjectedImagePlus = new ImagePlus("Z-Projected " + mOriginalImagePlus.getTitle(), vZProjectedStack);
//		vZProjectedImage.show();
    }

    private void filterTheInitialization(ImageStack aImageStack, int aInitPFIterations, int aFrameOfInit) {
        //
        // Preprocess the data several times on the first processor
        //
        Vector<Vector<float[]>> vInitStatesMemory = new Vector<Vector<float[]>>();
        vInitStatesMemory.add(copyStateVector(mStateVectors));

        // save a copy of the original sigma to restore it afterwards
        float[] vSigmaOfRWSave = new float[mSigmaOfRandomWalk.length];
        for (int vI = 0; vI < mSigmaOfRandomWalk.length; vI++) {
            vSigmaOfRWSave[vI] = mSigmaOfRandomWalk[vI];
        }

        for (int vR = 0; vR < aInitPFIterations; vR++) {
            IJ.showProgress(vR, aInitPFIterations);
            IJ.showStatus("Init progress " + vR / aInitPFIterations * 100 + "%");

            scaleSigmaOfRW(1f / (float) Math.pow(3, vR));
            DrawParticlesWithRW(mParticles);

            updateParticleWeights(aImageStack, aFrameOfInit);

            estimateStateVectors(mStateVectors, mParticles);

            resample(mParticles);

            vInitStatesMemory.add(copyStateVector(mStateVectors));
        }

        //restore the sigma vector
        for (int vI = 0; vI < mSigmaOfRandomWalk.length; vI++) {
            mSigmaOfRandomWalk[vI] = vSigmaOfRWSave[vI];
        }

//		if(true)
//			printStatesToWindow("Initialization states",vInitStatesMemory);

    }

    private void runParticleFilter(ImagePlus aImagePlus) {
        //
        // Copy the particles of the first frame if necessary(they were created while initialization of the filter)
        //
        if (mDoMonitorParticles) {
            mParticleMonitor.add(copyParticleVector(mParticles));
        }

         for (int vFrameIndex = mFrameOfInitialization; vFrameIndex <= mTrackTillFrameNb; vFrameIndex++) {
            ImageStack vCurrentFloatFrame = getAFrameCopy(aImagePlus, vFrameIndex);
            IJ.showProgress(vFrameIndex, aImagePlus.getNFrames());
            IJ.showStatus("Particle Filter in progress at frame: " + vFrameIndex);

            // save a copy of the original sigma to restore it afterwards
            float[] vSigmaOfRWSave = new float[mSigmaOfRandomWalk.length];
            for (int vI = 0; vI < mSigmaOfRandomWalk.length; vI++) {
                vSigmaOfRWSave[vI] = mSigmaOfRandomWalk[vI];
            }

            for (int vRepStep = 0; vRepStep < mRepSteps; vRepStep++) {
                if (vRepStep == 0) {
                    drawNewParticles(mParticles); //draw the particles at the appropriate position.
                } else {
                    scaleSigmaOfRW(1f / (float) Math.pow(3, vRepStep));//(1f - (float)vRepStep / (float)mRepSteps);
                    DrawParticlesWithRW(mParticles);
                }

                updateParticleWeights(vCurrentFloatFrame, vFrameIndex);

                estimateStateVectors(mStateVectors, mParticles);

                if (mDoResampling) {
                    if (!resample(mParticles)) {//further iterations are not necessary.
//						System.out.println("number of iterations needed at this frame: " + vRepStep);
                        break;
                    }
                }
                if (!mDoPrecisionOptimization) {
                    break; //do not repeat the filter on the first frame
                }
            }

            if (mDoMonitorParticles) {
                mParticleMonitor.add(copyParticleVector(mParticles));
            }

            //restore the sigma vector
            for (int vI = 0; vI < mSigmaOfRandomWalk.length; vI++) {
                mSigmaOfRandomWalk[vI] = vSigmaOfRWSave[vI];
            }
//			if(IJ.escapePressed()) { //TODO: doesnt work!
//				break;
//			}
            //save the new states
            mStateVectorsMemory.set(vFrameIndex - 1, copyStateVector(mStateVectors));
            //update the view

        }
        IJ.freeMemory();
    }

    /**
     * If one uses multiple iterations on a frame, this method scales the search
     * radius 'sigma'.
     *
     * @param vScaler
     */
    protected void scaleSigmaOfRW(float vScaler) {
        for (int vI = 0; vI < mSigmaOfRandomWalk.length; vI++) {
            mSigmaOfRandomWalk[vI] *= vScaler;
        }
    }

    private void createParticles(Vector<float[]> aStateVectors, Vector<Vector<float[]>> aParticles) {
        int vStateCounter = 0;
        if (aStateVectors.isEmpty()) {
            throw new IllegalArgumentException();
        }
        aParticles.clear();
        for (float[] vState : aStateVectors) {
            Vector<float[]> vParticleVector = new Vector<float[]>(mNbParticles);
//			float vPxWidthInNm = getPixelWidthInNm();
//			float vPxDepthInNm = getPixelDepthInNm();
            for (int vIndex = 0; vIndex < mNbParticles; vIndex++) {
                float[] vProposal = new float[vState.length + 1];
                for (int vI = 0; vI < vState.length; vI++) {
                    vProposal[vI] = vState[vI];
                }
//				Init the weight as a last dimension
                vProposal[vState.length] = 1f; //not 0!

                //Draw a the new proposal
//				drawFromProposalDistribution(vProposal, vPxWidthInNm, vPxDepthInNm);				

                //add the new particle
                vParticleVector.add(vProposal);
            }
            aParticles.add(vParticleVector);

            vStateCounter++;
        }
    }

    private void initMonitoring() {
        if (mDoMonitorIdealImage) {
            mIdealImageMonitorStack = new ImageStack(mWidth, mHeight);
        }
    }

    private void visualizeMonitors() {
        if (mDoMonitorIdealImage) {
            ImagePlus vIdealImagePlus = new ImagePlus("Summed up ideal image", mIdealImageMonitorStack);
            new StackWindow(vIdealImagePlus);
        }

        if (mDoMonitorParticles) {
            ImageStack vStack = new ImageStack(mWidth, mHeight);
            for (int vF = 0; vF < mOriginalImagePlus.getNFrames() * mRepSteps + 1; vF++) {
                vStack.addSlice("", new ByteProcessor(mWidth, mHeight));
            }
            ImagePlus vParticlesImage = new ImagePlus("Particles", vStack);
            new StackWindow(vParticlesImage, new ParticleMonitorCanvas(vParticlesImage));
        }
    }

    /**
     * Draws new particles for all the objects
     *
     */
    private void drawNewParticles(Vector<Vector<float[]>> aParticlesToRedraw) {
        //invoke this method here to not repeat it for every particle, pass it by argument
        //TODO: Better would probably be to scale the sigma vector from beginning; this would then introduce the problem that 
        //		the sampling of the particles cannot be treated for each dimension independently, which introduces
        //		an error.
        float vPxW = getPixelWidthInNm();
        float vPxD = getPixelDepthInNm();
        for (Vector<float[]> vObjectParticles : aParticlesToRedraw) {
            for (float[] vParticle : vObjectParticles) {
                drawFromProposalDistribution(vParticle, vPxW, vPxD);
            }
        }
    }

    /**
     *
     * @param aParticles set of parameters to resample.
     * @return true if resampling was performed, false if not.
     */
    private boolean resample(Vector<Vector<float[]>> aParticles) {
        for (Vector<float[]> vObjectParticles : aParticles) {
            //
            // First check if the threshold is smaller than Neff
            //
            int vDimOfState = aParticles.elementAt(0).elementAt(0).length - 1;
            float vNeff = 0;
            for (float[] vParticle : vObjectParticles) {
                vNeff += vParticle[vDimOfState] * vParticle[vDimOfState];
            }
            vNeff = 1 / vNeff;

            if (vNeff > mResamplingThreshold) {
//				System.out.println("no resampling");
                return false; //we won't do the resampling
            }
            //
            // Begin resampling
            //
//			System.out.println("Resampling");
            float VNBPARTICLES_1 = 1f / (float) mNbParticles;
            double[] vC = new double[mNbParticles + 1];
            vC[0] = 0;
            for (int vInd = 1; vInd <= mNbParticles; vInd++) {
                vC[vInd] = vC[vInd - 1] + vObjectParticles.elementAt(vInd - 1)[vDimOfState];
            }

            double vU = mRandomGenerator.nextFloat() * VNBPARTICLES_1;

            Vector<float[]> vFPParticlesCopy = copyStateVector(vObjectParticles);
            int vI = 0;
            for (int vParticleCounter = 0; vParticleCounter < mNbParticles; vParticleCounter++) {
                while (vU > vC[vI]) {
                    if (vI < mNbParticles) //this can happen due to numerical reasons
                    {
                        vI++;
                    }
                }
                for (int vK = 0; vK < vDimOfState; vK++) {
                    vObjectParticles.elementAt(vParticleCounter)[vK] = vFPParticlesCopy.elementAt(vI - 1)[vK];
                }
                vObjectParticles.elementAt(vParticleCounter)[vDimOfState] = VNBPARTICLES_1;
                vU += VNBPARTICLES_1;
            }
        }
        return true;
    }

    /**
     * Updates and normalizes the weights of the particles for all the objects
     *
     * @param aObservationImage The picture of the next time step
     * @param aFrameIndex: The frame index corresponding to picture in
     * aObservationStack
     */
    private void updateParticleWeights(ImageStack aObservationStack, int aFrameIndex) {
        int vObjectIndex = 0;
        for (Vector<float[]> vObjectParticles : mParticles) {
            vObjectIndex++;
            float vSumOfWeights = 0;
            int vDimOfState = vObjectParticles.elementAt(0).length - 1;

            //
            // Calculate the likelihoods for each particle and save the biggest one
            //
            float[] vLogLikelihoods = new float[mNbParticles];
            float vMaxLogLikelihood = Float.NEGATIVE_INFINITY;
            boolean[][][] vBitmap = generateParticlesIntensityBitmap_3D(vObjectParticles, mWidth, mHeight, mNSlices);
            Thread[] vThreads = new Thread[mNbThreads];
            for (int vT = 0; vT < mNbThreads; vT++) {
                vThreads[vT] = new ParallelizedLikelihoodCalculator(aObservationStack, vBitmap, aFrameIndex, vLogLikelihoods, vObjectParticles);
            }
            for (int vT = 0; vT < mNbThreads; vT++) {
                vThreads[vT].start();
            }
            //wait for the threads to end.
            for (int vT = 0; vT < mNbThreads; vT++) {
                try {
                    vThreads[vT].join();
                } catch (InterruptedException aIE) {
                    IJ.showMessage("Not all particles calculated, the tracking might be wrong.");
                }
            }

            for (int vI = 0; vI < vObjectParticles.size(); vI++) {
                if (vLogLikelihoods[vI] > vMaxLogLikelihood) {
                    vMaxLogLikelihood = vLogLikelihoods[vI];
                }
            }
            mMaxLogLikelihood[aFrameIndex - 1] = vMaxLogLikelihood;
            //
            // Iterate again and update the weights
            //
            int vI = 0;
            for (float[] vParticle : vObjectParticles) {
                vLogLikelihoods[vI] -= vMaxLogLikelihood;
                vParticle[vDimOfState] = vParticle[vDimOfState] * (float) Math.exp(vLogLikelihoods[vI]);
                vSumOfWeights += vParticle[vDimOfState];
                vI++;
            }
            //
            // Iterate again and normalize the weights
            //
            if (vSumOfWeights == 0.0f) { //can happen if the winning particle before had a weight of 0.0
                for (float[] vParticle : vObjectParticles) {
                    vParticle[vDimOfState] = 1.0f / (float) mNbParticles;
                }
            } else {
                for (float[] vParticle : vObjectParticles) {
                    vParticle[vDimOfState] /= vSumOfWeights;
                }
            }
        }
    }

    /**
     * Estimates all state vectors from the particles and their weights
     */
    private void estimateStateVectors(Vector<float[]> aStateVectors, Vector<Vector<float[]>> aParticles) {
        for (int vFPIndex = 0; vFPIndex < aStateVectors.size(); vFPIndex++) {
            float[] vState = aStateVectors.get(vFPIndex);
            Vector<float[]> vObjectParticles = aParticles.get(vFPIndex);

            /*
             * Set the old state to 0
             */
            for (int vI = 0; vI < vState.length; vI++) {
                vState[vI] = 0f;
            }

            for (float[] vParticle : vObjectParticles) {
                for (int vDim = 0; vDim < vState.length; vDim++) {
                    vState[vDim] += vParticle[vState.length] * vParticle[vDim];
                }
            }
        }
    }
    protected boolean[][][] mIntensityBitmap = null;

    /**
     * Override this method to speed up the algorithm (useful in 3D images). It
     * may be useful to fill the member
     * <code>mIntensityBitmap</code> in the
     * <code>generateIdealImage_3D(...)</code> method and return it here.
     *
     * @param aSetOfParticles
     * @return
     */
    protected boolean[][][] generateParticlesIntensityBitmap_3D(Vector<float[]> aSetOfParticles, int aW, int aH, int aS) {
        if (mIntensityBitmap == null) {
            mIntensityBitmap = new boolean[aS][aH][aW];
            for (int vZ = 0; vZ < aS; vZ++) {
                for (int vY = 0; vY < aH; vY++) {
                    for (int vX = 0; vX < aW; vX++) {
                        mIntensityBitmap[vZ][vY][vX] = true;
                    }
                }
            }
        }
        return mIntensityBitmap;
    }

    /**
     * Generates an artificial 3D image. Override this method. Do not forget to
     * add background(minimal value of a voxel is 1.0).
     *
     * @param aW: The width of the image to generate
     * @param aH: The height of the image to generate
     * @param aS: The number of slices of the image to generate
     * @param aParticle: the particle that describes the state
     * @param aBackground: background intensity to add
     * @param aPxDepthInNm
     * @param aPxWidthInNm
     * @return a 3D image
     */
    protected abstract float[][][] generateIdealImage_3D(int aW, int aH, int aS, float[] aParticle, int aBackground, float aPxWidthInNm, float aPxDepthInNm);

    /**
     * The plugin can draw on the canvas.
     *
     * @param aG: The graphics object to draw on.
     * @param aMagnification: The magnification factor used by the user.
     * @param aActiveFrame: The frame currently selected by the user.
     */
    protected abstract void paintOnCanvas(Graphics aG, double aMagnification, int aActiveFrame);

    /**
     * Calculates the likelihood by multipling the poissons marginals around a
     * particle given a image(optimal image)
     *
     * @param aImagePlus: The observed image
     * @param aFrame: The frame index 1<=n<=NSlices(to read out the correct
     * substack from aStack
     * @param aGivenImage: a intensity array, the 'measurement'
     * @param aBitmap: Pixels which are not set to true in the bitmap are pulled
     * out(proportionality)
     * @return the likelihood for the image given
     */
    private float calculateLogLikelihood_3D(float[][] aStackProcs, int aFrame, float[][][] aGivenImage, boolean[][][] aBitmap) //ImageStack aImageStack){
    {
        float vLogLikelihood = 0;
        //we need all processors anyway. Profiling showed that the method getProcessor needs a lot of time. Store them
        //in an Array.

//		long vTime1 = System.currentTimeMillis();
        for (int vZ = 0; vZ < mNSlices; vZ++) {
            for (int vY = 0; vY < mHeight; vY++) {
                for (int vX = 0; vX < mWidth; vX++) {
                    if (aBitmap[vZ][vY][vX]) {
                        vLogLikelihood += -aGivenImage[vZ][vY][vX] + (float) aStackProcs[vZ][vY * mWidth + vX] * (float) Math.log(aGivenImage[vZ][vY][vX]);
                        if (Float.isNaN(vLogLikelihood)) {
                            System.out.println("NAN at vz = " + vZ + ", vY = " + vY + ", vX = " + vX);
                        }
                    }
                }
            }
        }
//		System.out.println("used time for loglik = " + (System.currentTimeMillis() - vTime1));
        //IJ.showStatus("likelihood finshed");
        return vLogLikelihood;
    }

    /**
     * Override this method to include the apriori knowledge of the dynamics of
     * your system.
     *
     * @param aParticle A possible state.
     * @param aPxWidthInNm Pixel width in nano meter to not use the getter
     * function since it is called very often.
     * @param aPxDepthInNm Pixel depth in nano meter to not use the getter
     * function since it is called very often.
     */
    abstract protected void drawFromProposalDistribution(float[] aParticle, float aPxWidthInNm, float aPxDepthInNm);

    /**
     * Converts a slice index to a frame index using the parameters of the image
     * defined by the user.
     *
     * @param aSlice
     * @return frame index.
     */
    protected int sliceToFrame(int aSlice) {
        if (aSlice < 1) {
            System.err.println("wrong argument in particle filter in SliceToFrame: < 1");
        }
        return (int) (aSlice - 1) / mOriginalImagePlus.getNSlices() + 1;
    }

    private void DrawParticlesWithRW(Vector<Vector<float[]>> aParticles) {
        float vPxW = getPixelWidthInNm();
        float vPxD = getPixelDepthInNm();
        for (Vector<float[]> vFPParticles : aParticles) {
            for (float[] vP : vFPParticles) {
                randomWalkProposal(vP, vPxW, vPxD);
            }
        }
    }

    /**
     * Draws new particles, the last entry of the argument is supposed to be the
     * weight and remains unchanged.
     *
     * @param aParticle A Array with state vector entries + a weight
     * entry(remains unchanged)
     */
    private void randomWalkProposal(float[] aParticle, float aPxWidthInNm, float aPxDepthInNm) {
        for (int aI = 0; aI < aParticle.length - 1; aI++) {
            aParticle[aI] += (float) mRandomGenerator.nextGaussian() * mSigmaOfRandomWalk[aI];
        }
    }

    /**
     * Checks if the movie/image is saved and returns the possible result file.
     *
     * @return null if the file is not saved, else it returns the file(even if
     * it does not exist).
     */
    protected File getResultFile() {
        FileInfo vFI = mOriginalImagePlus.getOriginalFileInfo();
        if (vFI == null) {
            return null;
        }
        String vResFileName = new String(vFI.fileName);
        int vLastInd = vResFileName.lastIndexOf(".");
        if (vLastInd != -1) {
            vResFileName = vResFileName.substring(0, vLastInd);
        }
        vResFileName = vResFileName.concat(RESULT_FILE_SUFFIX);

        File vResFile = new File(vFI.directory, vResFileName);
        return vResFile;
    }

    /**
     * @return null if the movie/image is not saved on disk. Otherwise the init
     * file is returned.
     * @see getResultFile()
     */
    protected File getInitFile() {
        FileInfo vFI = mOriginalImagePlus.getOriginalFileInfo();
        if (vFI == null) {
            return null;
        }
        String vFileName = new String(vFI.fileName);
        int vLastInd = vFileName.lastIndexOf(".");
        vLastInd = vFileName.lastIndexOf(".");
        if (vLastInd != -1) {
            vFileName = vFileName.substring(0, vLastInd);
        }
        vFileName = vFileName.concat(INIT_FILE_SUFFIX);

        File vInitFile = new File(vFI.directory, vFileName);
        return vInitFile;

    }

    /**
     * Writes the initialization to disk. Using this information, batch
     * processing movies can be done using macros.
     *
     * @param aFile
     * @return true if successful, false if not.
     */
    protected boolean writeInitFile(File aFile) {
        BufferedWriter vW = null;
        try {
            vW = new BufferedWriter(new FileWriter(aFile));
            for (float[] vState : mStateVectors) {
                String vS = mFrameOfInitialization + " ";
                for (int vI = 0; vI < vState.length; vI++) {
                    vS += vState[vI] + " ";
                }
                vW.write(vS + "\n");
            }
        } catch (IOException aIOE) {
            aIOE.printStackTrace();
            return false;
        } finally {
            try {
                vW.close();
            } catch (IOException aIOE) {
                aIOE.printStackTrace();
                return false;
            }
        }
        return true;
    }

    /**
     * Writes the results to the result file.
     *
     * @see getResultFile()
     * @param aFile
     * @return true if successful, false if not.
     */
    protected boolean writeResultFile(File aFile) {
        BufferedWriter vW = null;
        try {
            vW = new BufferedWriter(new FileWriter(aFile));
            vW.write(generateOutputString(mStateVectorsMemory, ",", true));
        } catch (IOException aIOE) {
            aIOE.printStackTrace();
            return false;
        } finally {
            try {
                vW.close();
            } catch (IOException aIOE) {
                aIOE.printStackTrace();
                return false;
            }
        }
        return true;
    }

    /**
     * Reads the init file. The parameters read are then stored in the
     * statevector member.
     *
     * @param aFile
     * @see getInitFile()
     * @return true if successful. false if not.
     */
    protected boolean readInitFile(File aFile) {
        BufferedReader vR = null;
        try {
            vR = new BufferedReader(new FileReader(aFile));
        } catch (FileNotFoundException aFNFE) {
            return false;
        }
        String vLine;
        try {
            while ((vLine = vR.readLine()) != null) {
                if (vLine.startsWith("#")) {
                    continue; //comment
                }
                if (vLine.matches("(\\s)*")) {
                    continue; //empty line
                }
                Pattern vPattern = Pattern.compile("(,|\\||;|\\s)");
                String[] vPieces = vPattern.split(vLine);
//				if (vPieces.length < mStateVectors.firstElement().length) continue; //not the right line
                try {
                    mFrameOfInitialization = Integer.parseInt(vPieces[0]);
                } catch (NumberFormatException aNFE) {
                    continue;
                }
                float[] vState = new float[vPieces.length - 1];
                for (int i = 1; i < vPieces.length; i++) {
                    float vValue = 0f;
                    try {
                        vValue = Float.parseFloat(vPieces[i]);
                        vState[i - 1] = vValue;
                    } catch (NumberFormatException aNFE) {
                        continue; //perhaps there is another matching line
                    }
                }
                mStateVectors.add(vState);

            }
        } catch (IOException aIOE) {
            aIOE.printStackTrace();
            return false;
        } finally {
            try {
                vR.close();
            } catch (IOException aIOE) {
                aIOE.printStackTrace();
                return false;
            }
        }

        return true;
    }

    /**
     * Reads the result file. The parameters read are then stored in the
     * statevectormemory member.
     *
     * @param aFile
     * @see getInitFile()
     * @return true if successful. false if not.
     */
    protected boolean readResultFile(File aFile) {
        BufferedReader vR = null;
        try {
            vR = new BufferedReader(new FileReader(aFile));
        } catch (FileNotFoundException aFNFE) {
            return false;
        }
        String vLine;
        try {
            while ((vLine = vR.readLine()) != null) {
                if (vLine.startsWith("#")) {
                    continue; //ignore
                }
                if (vLine.matches("(\\s)*")) {
                    continue; //empty line
                }
                Pattern vPattern = Pattern.compile("(,|\\||;|\\s)");
                String[] vPieces = vPattern.split(vLine);

                if (vPieces[0].equalsIgnoreCase("frame")) {
                    mDimensionsDescription = new String[vPieces.length - 1];
                    for (int vP = 1; vP < vPieces.length; vP++) {
                        mDimensionsDescription[vP - 1] = vPieces[vP];
                    }
                }
                //if (vPieces.length < mStateVectors.firstElement().length) continue; //not the right line
                int vFrame = 0;
                try {
                    vFrame = Integer.parseInt(vPieces[0]) - 1;
                    if (vFrame < 0) {
                        IJ.showMessage("Warning", "Unproper result file");
                        continue;
                    }
                } catch (NumberFormatException aNFE) {
                    continue;
                }
                Vector<float[]> vFrameStates;
                if (mStateVectorsMemory.elementAt(vFrame) == null) {
                    vFrameStates = new Vector<float[]>();
                } else {
                    vFrameStates = mStateVectorsMemory.elementAt(vFrame);
                }
                float[] vFrameState = new float[vPieces.length - 1];

                for (int vP = 1; vP < vPieces.length; vP++) {
                    float vValue = 0f;
                    try {
                        vValue = Float.parseFloat(vPieces[vP]);
                        vFrameState[vP - 1] = vValue;
                    } catch (NumberFormatException aNFE) {
                        continue;//perhaps there is another matching line
                    }
                }
                vFrameStates.add(vFrameState);
                mStateVectorsMemory.setElementAt(vFrameStates, vFrame);

            }
        } catch (IOException aIOE) {
            aIOE.printStackTrace();
            return false;
        } finally {
            try {
                vR.close();
            } catch (IOException aIOE) {
                aIOE.printStackTrace();
                return false;
            }
        }

        return true;
    }

    /**
     * This method may be overrided.
     *
     * @param aStateVectorsMem
     * @return the string to print in a file or in a window(if the movie is not
     * saved).
     */
    protected String generateOutputString(Vector<Vector<float[]>> aStateVectorsMem, String aDelimiter, boolean aPrintHeader) {
        //
        // Get the number of objects.o
        //
        int vNObjects = 0;
        for (int vI = 0; vI < aStateVectorsMem.size(); vI++) {
            if (aStateVectorsMem.elementAt(vI) != null) {
                vNObjects = aStateVectorsMem.elementAt(vI).size();
                break;
            }
        }
        String vOut = "";

        if (aPrintHeader) {
            vOut += "frame" + aDelimiter;
            for (String vDim : mDimensionsDescription) {
                vOut += vDim + aDelimiter;
            }
        }

        for (int vI = 0; vI < vNObjects; vI++) {
            int vFrameC = 0;
            for (Vector<float[]> vFrame : aStateVectorsMem) {
                if (vFrame == null) {
                    vFrameC++;
                    continue;
                }
                vFrameC++;
                vOut += "\n" + vFrameC + aDelimiter;

                for (float vV : vFrame.elementAt(vI)) {
                    vOut += vV + aDelimiter;
                }

//				vOut += mMaxLogLikelihood[vFrameC-1] + ",";				
            }
            vOut += "\n";
        }
        return vOut;
    }

    /**
     * Generates a text window to display all the information in
     * <code>aStateVectorsMem</code>
     *
     * @param aTitle Title of the window
     * @param aStateVectorsMem
     */
    protected void printStatesToWindow(String aTitle, Vector<Vector<float[]>> aStateVectorsMem) {
        String vData = generateOutputString(aStateVectorsMem, "\t", true);
        int vFirstLineEndIndex = vData.indexOf("\n");
        if (vFirstLineEndIndex < 0) {
            IJ.showMessage("Empty.");
            return;
        }
        new TextWindow(aTitle, vData.substring(0, vFirstLineEndIndex),
                vData.substring(vFirstLineEndIndex + 1, vData.length()), 400, 400);
    }

    /**
     * Reads in the parameters used by the base class.
     *
     * @return True, if cancelled false.
     */
    protected boolean getUserDefinedParams() {
        mDoPrecisionOptimization = getMDoPrecisionOptimization();

        GenericDialog vGenericDialog = new GenericDialog("Particle filtering parameters", IJ.getInstance());
        vGenericDialog.addNumericField("# of particles ~ quality", mNbParticles, 0);
        if (mDoPrecisionOptimization) {
            vGenericDialog.addNumericField("Maximal iterations on a frame", mRepSteps, 0);
        }
        vGenericDialog.addNumericField("Number of iterations 1st frame", mInitRWIterations, 0);
        vGenericDialog.addNumericField("Background Intensity", mBackground, 0);
        vGenericDialog.addNumericField("Wavelength in nm", mWavelengthInNm, 0);
        vGenericDialog.addNumericField("Numerical apparture", mNA, 2);
        vGenericDialog.addNumericField("refractive index(medium)", mn, 2);
        vGenericDialog.addNumericField("Track till frame", mTrackTillFrameNb, 0);

        vGenericDialog.showDialog();

        if (vGenericDialog.wasCanceled()) {
            return false;
        }

        mNbParticles = (int) vGenericDialog.getNextNumber();
        if (mDoPrecisionOptimization) {
            mRepSteps = (int) vGenericDialog.getNextNumber();
        }
        mInitRWIterations = (int) vGenericDialog.getNextNumber();
        mBackground = (int) vGenericDialog.getNextNumber(); //not necessary
        setMBackground(mBackground);
        mWavelengthInNm = (int) vGenericDialog.getNextNumber();
        mNA = (float) vGenericDialog.getNextNumber();
        mn = (float) vGenericDialog.getNextNumber();
        mTrackTillFrameNb = (int) vGenericDialog.getNextNumber();

        return true;
    }

    /**
     * Shows the dialog to enter the specific parameters.
     *
     * @return false if cancelled, else true.
     */
    protected boolean showParameterDialog() {
        return true;
    }

    /**
     * Copies a
     * <code>Vector&lt;float[]&gt;</code> data structure. Used to copy the state
     * vector here.
     *
     * @param aOrig
     * @return the copy.
     */
    public static Vector<float[]> copyStateVector(Vector<float[]> aOrig) {
        Vector<float[]> vResVector = new Vector<float[]>(aOrig.size());
        for (float[] vA : aOrig) {
            float[] vResA = new float[vA.length];
            for (int vI = 0; vI < vA.length; vI++) {
                vResA[vI] = vA[vI];
            }
            vResVector.add(vResA);
        }
        return vResVector;
    }

    /**
     * Copies a
     * <code>Vector&lt;Vector&lt;float[]&gt;&gt;</code> data structure. Used to
     * copy the particle vector here.
     *
     * @param aOrig
     * @return the copy.
     */
    public static Vector<Vector<float[]>> copyParticleVector(Vector<Vector<float[]>> aOrig) {
        Vector<Vector<float[]>> vResVector = new Vector<Vector<float[]>>(aOrig.size());
        for (Vector<float[]> vP : aOrig) {
            vResVector.add(copyStateVector(vP));
        }
        return vResVector;
    }

    /**
     * Recursively searches the brightest voxel in the neighborhood. Might be
     * used for the initialization.
     *
     * @param aStartX
     * @param aStartY
     * @param aStartZ
     * @param aImageStack
     * @return a int array with 3 entries: x,y and z coordinate.
     */
    public static int[] searchLocalMaximumIntensityWithSteepestAscent(int aStartX, int aStartY, int aStartZ, ImageStack aImageStack) {
        int[] vRes = new int[]{aStartX, aStartY, aStartZ};
        float vMaxValue = aImageStack.getProcessor(aStartZ).getPixelValue(aStartX, aStartY);
        for (int vZi = -1; vZi < 2; vZi++) {
            if (aStartZ + vZi > 0 && aStartZ + vZi <= aImageStack.getSize()) {
                for (int vXi = -1; vXi < 2; vXi++) {
                    for (int vYi = -1; vYi < 2; vYi++) {
                        if (aImageStack.getProcessor(aStartZ + vZi).getPixelValue(aStartX + vXi, aStartY + vYi) > vMaxValue) {
                            vMaxValue = aImageStack.getProcessor(aStartZ + vZi).getPixelValue(aStartX + vXi, aStartY + vYi);
                            vRes[0] = aStartX + vXi;
                            vRes[1] = aStartY + vYi;
                            vRes[2] = aStartZ + vZi;
                        }
                    }
                }
            }
        }
        if (vMaxValue > aImageStack.getProcessor(aStartZ).getPixelValue(aStartX, aStartY)) {
            return searchLocalMaximumIntensityWithSteepestAscent(vRes[0], vRes[1], vRes[2], aImageStack);
        }
        return vRes;
    }

    /**
     * Add the intensities of 2 2D arrays.
     *
     * @param aResult here the first image is stored in.
     * @param aImageToAdd a Image that is added to
     * <code>aResult</code>
     */
    public static void addImage(float[][] aResult, float[][] aImageToAdd) {
        int vIMax = Math.min(aResult.length, aImageToAdd.length);
        int vJMax = Math.min(aResult[0].length, aImageToAdd[0].length);
        for (int vI = 0; vI < vIMax; vI++) {
            for (int vJ = 0; vJ < vJMax; vJ++) {
                aResult[vI][vJ] += aImageToAdd[vI][vJ];
            }
        }
    }

    /**
     * Sets all values in the array to
     * <code>aValue</code>
     *
     * @param aArray
     * @param aValue
     */
    public static void initArrayToValue(float[][] aArray, float aValue) {
        for (int vI = 0; vI < aArray.length; vI++) {
            for (int vJ = 0; vJ < aArray[0].length; vJ++) {
                aArray[vI][vJ] = aValue;
            }
        }
    }

    /**
     * Returns a copy of a single frame. Note that the properties of the
     * ImagePlus have to be correct
     *
     * @param aMovie
     * @param aFrameNumber
     * @return The frame copy.
     */
    public static ImageStack getAFrameCopy(ImagePlus aMovie, int aFrameNumber) {
        if (aFrameNumber > aMovie.getNFrames() || aFrameNumber < 1) {
            throw new IllegalArgumentException();
        }
        int vS = aMovie.getNSlices();
        return getSubStackFloatCopy(aMovie.getStack(), (aFrameNumber - 1) * vS + 1, aFrameNumber * vS);
    }

    /**
     * Rerurns a copy of a substack (i.e.frames)
     *
     * @param aImageStack: the stack to crop
     * @param aStartPos: 1 &le; aStartPos &le; aImageStack.size()
     * @param aEndPos: 1 &le; aStartPos &le; aEndPos &le; aImageStack.size()
     * @return a Copy of the supstack
     */
    public static ImageStack getSubStackFloatCopy(ImageStack aImageStack, int aStartPos, int aEndPos) {
        ImageStack res = new ImageStack(aImageStack.getWidth(), aImageStack.getHeight());
        if (!(aStartPos < 1 || aEndPos < 0)) {
            for (int vI = aStartPos; vI <= aEndPos; vI++) {
                res.addSlice(aImageStack.getSliceLabel(vI), aImageStack.getProcessor(vI).convertToFloat().duplicate());
            }
        }
        return res;
    }

    /**
     *
     * @param aImageStack: the stack to crop
     * @param aStartPos: 1 &le; aStartPos &le; aImageStack.size()
     * @param aEndPos: 1 &le; aStartPos &le; aEndPos &le; aImageStack.size()
     * @return
     */
    public static ImageStack getSubStackFloat(ImageStack aImageStack, int aStartPos, int aEndPos) {
        ImageStack res = new ImageStack(aImageStack.getWidth(), aImageStack.getHeight());
        if (!(aStartPos < 1 || aEndPos < 0)) {
            for (int vI = aStartPos; vI <= aEndPos; vI++) {
                res.addSlice(aImageStack.getSliceLabel(vI), aImageStack.getProcessor(vI).convertToFloat());
            }
        }
        return res;
    }
    private int mControllingParticleIndex = 0;

    @SuppressWarnings("serial")
    private class DrawCanvas extends ImageCanvas {

        public DrawCanvas(ImagePlus aImagePlus) {
            super(aImagePlus);
        }

        public void paint(Graphics aG) {
            super.paint(aG);
            paintOnCanvas(aG, magnification, mZProjectedImagePlus.getCurrentSlice());
        }
    }

    @SuppressWarnings("serial")
    private class ParticleMonitorCanvas extends ImageCanvas {

        ImagePlus mImagePlus;

        public ParticleMonitorCanvas(ImagePlus aImagePlus) {
            super(aImagePlus);
            mImagePlus = aImagePlus;
        }

        public void paint(Graphics aG) {
            super.paint(aG);
            try {
                if (mParticleMonitor.size() >= mImagePlus.getCurrentSlice()) {
                    for (Vector<float[]> vObjectParticle : mParticleMonitor.elementAt(mImagePlus.getCurrentSlice() - 1)) {
                        for (float[] vParticle : vObjectParticle) {
                            paintParticleOnCanvas(aG, vParticle, magnification);
                        }
                    }
                }
            } catch (java.lang.NullPointerException vE) {
                //do nothing
            }
        }
    }

    /**
     * To override if you would like to visualize the particles in a separate
     * window.
     *
     * @param aParticle
     * @param aMagnification
     */
    protected void paintParticleOnCanvas(Graphics aG, float[] aParticle, double aMagnification) {
        return;
    }

    private class ParallelizedLikelihoodCalculator extends Thread {
        //
        // Most of this members are used to speed up the algorithm, so that they have not to be
        // evaluated for each particle but only for each thread.
        //

        float[] mResultArray;//we only write to the array, there should be no conflicts
        Vector<float[]> mParticles;
        ImageStack mObservedImage;
        float[][] mStackProcs;
        boolean[][][] mBitmap;
        int mFrameIndex;
        float mPxWidthInNm, mPxDepthInNm;

        /**
         * Calculates likelihoods for particles given a image and writes them in
         * the result array. There are two options: 1. LeastSquares
         * likelihoods(negative!) and 2. Poisson LOG(!) Likelihoods. DO FIRST
         * CONSTRUCT ALL THREADS BEFORE RUNNING THE FIRST ONE!
         *
         * @param aImageStack: The frame to operate on, i.e. the observed image
         * @param aRusultArray: the results are written in this array in the
         * same order the particles are in
         * <code>aParticles</code>
         * @param aParticles: the particles to score
         *
         */
        public ParallelizedLikelihoodCalculator(ImageStack aImageStack, boolean[][][] aBitmap, int aFrameIndex, float[] aRusultArray, Vector<float[]> aParticles) {
            mResultArray = aRusultArray;
            mParticles = aParticles;
            mObservedImage = aImageStack;
            mBitmap = aBitmap;
            mFrameIndex = aFrameIndex;
            // The next line is only ok if the threads are first constructed and run afterwards!
            mControllingParticleIndex = 0;
            //this members speeds the algorithm drastically up since we only have to 
            //invoke getProcessor() once per thread instead of for each particle.
            mStackProcs = new float[mNSlices][];
            mPxWidthInNm = getPixelWidthInNm();
            mPxDepthInNm = getPixelDepthInNm();
            for (int vZ = 0; vZ < mNSlices; vZ++) {
                mStackProcs[vZ] = (float[]) mObservedImage.getProcessor(vZ + 1).getPixels();
            }

        }

        public void run() {
            int vI;
            while ((vI = getNewParticleIndex()) != -1) {
//								if(vI == 2) {
//									System.out.println("run stop");
//								}
                //get the particle
                float[] vParticle = mParticles.elementAt(vI);
                //calculate ideal image
                float[][][] vIdealImage = generateIdealImage_3D(mWidth,
                        mHeight,
                        mNSlices,
                        vParticle,
                        (int) (mBackground + .5),
                        mPxWidthInNm,
                        mPxDepthInNm);

                //calculate likelihood
                mResultArray[vI] = calculateLogLikelihood_3D(mStackProcs, mFrameIndex, vIdealImage, mBitmap);
//				if(Float.isNaN(mResultArray[vI])){
//					System.out.println("NAN found! at particle index vI = " + vI);
//				}

            }
        }

        synchronized int getNewParticleIndex() {
            if (mControllingParticleIndex < mNbParticles && mControllingParticleIndex >= 0) {
                mControllingParticleIndex++;
                return mControllingParticleIndex - 1;
            }
            mControllingParticleIndex = -1;
            return -1;
        }
    }

    /**
     * Method to override optionally. The body of the method in the base class
     * is empty.
     *
     */
    protected void mousePressed(int aX, int aY) {
    }

    /**
     * Method to override optionally. The body of the method in the base class
     * is empty.
     *
     */
    protected void mouseClicked(int aX, int aY) {
    }

    /**
     * Method to override optionally. The body of the method in the base class
     * is empty.
     *
     * @param aEvent
     */
    protected void mouseEntered(MouseEvent aEvent) {
    }

    /**
     * Method to override optionally. The body of the method in the base class
     * is empty.
     *
     * @param aEvent
     */
    protected void mouseExited(MouseEvent aEvent) {
    }

    /**
     * Method to override optionally. The body of the method in the base class
     * is empty.
     *
     */
    protected void mouseReleased(int aX, int aY) {
    }

    /**
     * May be overrided optionally. The implementation of the base class runs
     * the algorithm.
     */
    protected void calcFromHereButtonPressed() {
        //relaunch the filter from the current slice on with the current initialization if there is one
        if (mStateOfFilter != STATE_OF_FILTER.VISUALIZING && mStateOfFilter != STATE_OF_FILTER.READY_TO_RUN) {
            IJ.showMessage("No initialization done yet.");
            return;
        }
        if (mStateOfFilter == STATE_OF_FILTER.VISUALIZING) { //restart from corrected frame
            if (mStateVectorsMemory.elementAt(mZProjectedImagePlus.getCurrentSlice() - 1) != null) {
                //first setup the state vector then run
                mStateVectors = copyStateVector(mStateVectorsMemory.elementAt(mZProjectedImagePlus.getCurrentSlice() - 1));
                mFrameOfInitialization = mZProjectedImagePlus.getCurrentSlice();
                PFTracking3D.this.run(new FloatProcessor(1, 1));
                mStateOfFilter = STATE_OF_FILTER.RUNNING;
            } else {
                IJ.showMessage("No initialization in this frame.");
            }
        }
        if (mStateOfFilter == STATE_OF_FILTER.READY_TO_RUN) {
            //we're sure that there is a correct initialization at a certain frame
            mStateOfFilter = STATE_OF_FILTER.RUNNING;
            PFTracking3D.this.run(new FloatProcessor(1, 1));
            mStateOfFilter = STATE_OF_FILTER.VISUALIZING;
        }
    }

    /**
     * May be overrided optionally. The implementation of the base class sets
     * the state of the plugin to
     * <code>INIT</code>.
     */
    protected void initializeWithMouseButtonPressed() {
        mStateOfFilter = STATE_OF_FILTER.INIT;
    }

    /**
     * May be overrided optionally. The implementation of the base class saves
     * the initialization values in a separate file. This enables batch
     * processing.
     */
    protected void saveInitButtonPressed() {
        if (mStateOfFilter != STATE_OF_FILTER.VISUALIZING && mStateOfFilter != STATE_OF_FILTER.READY_TO_RUN) {
            IJ.showMessage("No initialization done yet. Please initialize.");
            return;
        }
        if (mOriginalImagePlus.getOriginalFileInfo() == null) {
            IJ.showMessage("It seems that the movie was not saved. Save the movie in a directory with 'write' permission first.");
            return;
        }
        mFrameOfInitialization = mZProjectedImagePlus.getCurrentSlice();
        //write out the init positions in a file in the same directory

        writeInitFile(getInitFile());
    }

    /**
     * May be overrided optionally. The implementation of the base class deletes
     * all values in RAM and on the Disk(Init-value file and result file).
     * Afterwards, a new initialization has to be done.
     */
    protected void deleteAllButtonPressed() {
        //clear the memory vectors(also zmem) and repaint the window.
        mStateVectorsMemory.clear();
        mStateVectors.clear();
        mStateVectorsMemory.setSize(mNFrames);
        if (getInitFile() != null) {
            getInitFile().delete();
        }
        if (getResultFile() != null) {
            getResultFile().delete();
        }
        mStateOfFilter = STATE_OF_FILTER.WAITING;
        mZProjectedImagePlus.repaintWindow();
    }

    /**
     * May be overrided optionally. The implementation of the base class
     * rewrites the result file on the HD.
     */
    protected void saveCorrectionButtonPressed() {
        if (mOriginalImagePlus.getOriginalFileInfo() == null) {
            IJ.showMessage("It seems that the movie was not saved. Save the movie in a directory with 'write' permission first.");
            return;
        }
        //write out the positions in the res file
        writeResultFile(getResultFile());
    }

    @SuppressWarnings("serial")
    private class TrajectoryStackWindow extends StackWindow implements ActionListener, MouseListener {
//		private static final long serialVersionUID = 1L;

        private Button mCalcFromHereButton;
        private Button mMouseInitializationButton;
        private Button mSaveInitButton;
        private Button mDeleteAllButton;
        private Button mSaveCorrectionButton;
        private Button mChangeParametersButton;
        private Button mShowResultsButton;

        /**
         * Constructor. <br>Creates an instance of TrajectoryStackWindow from a
         * given
         * <code>ImagePlus</code> and
         * <code>ImageCanvas</code> and a creates GUI panel. <br>Adds this class
         * as a
         * <code>MouseListener</code> to the given
         * <code>ImageCanvas</code>
         *
         * @param aImagePlus
         * @param aImageCanvas
         */
        private TrajectoryStackWindow(ImagePlus aImagePlus, ImageCanvas aImageCanvas) {
            super(aImagePlus, aImageCanvas);
            aImageCanvas.addMouseListener(this);
            addPanel();
        }

        /**
         * Adds a Panel with filter options button in it to this window
         */
        private void addPanel() {
            Panel vButtonPanel = new Panel(new GridLayout(4, 3));
            mChangeParametersButton = new Button("Change parameters...");
            mCalcFromHereButton = new Button("Calculate!");
            mSaveInitButton = new Button("Save init position");
            mDeleteAllButton = new Button("Delete all");
            mSaveCorrectionButton = new Button("Write data to disk");
            mMouseInitializationButton = new Button("Initialize with mouse");
            mShowResultsButton = new Button("Show Results");
            mCalcFromHereButton.addActionListener(this);
            mSaveInitButton.addActionListener(this);
            mDeleteAllButton.addActionListener(this);
            mMouseInitializationButton.addActionListener(this);
            mSaveCorrectionButton.addActionListener(this);
            mChangeParametersButton.addActionListener(this);
            mShowResultsButton.addActionListener(this);

            vButtonPanel.add(mCalcFromHereButton);
            vButtonPanel.add(mSaveInitButton);
            vButtonPanel.add(mSaveCorrectionButton);
            vButtonPanel.add(mDeleteAllButton);
            vButtonPanel.add(mMouseInitializationButton);
            vButtonPanel.add(mSaveCorrectionButton);
            vButtonPanel.add(mChangeParametersButton);
            vButtonPanel.add(mShowResultsButton);
            add(vButtonPanel);
            pack();
            Dimension screen = Toolkit.getDefaultToolkit().getScreenSize();
            Point loc = getLocation();
            Dimension size = getSize();
            if (loc.y + size.height > screen.height) {
                getCanvas().zoomOut(0, 0);
            }
        }

        /**
         * Defines the action taken upon an
         * <code>ActionEvent</code> triggered from buttons that have class
         * <code>TrajectoryStackWindow</code> as their action listener: <br><code>Button filter_length</code>
         *
         * @see
         * java.awt.event.ActionListener#actionPerformed(java.awt.event.ActionEvent)
         */
        public synchronized void actionPerformed(ActionEvent aEvent) {
            Object vButton = aEvent.getSource();

            //
            // Show result window
            //
            if (vButton == mShowResultsButton) {
                PFTracking3D.this.printStatesToWindow("Results", mStateVectorsMemory);
            }

            //
            // Change param button
            //
            if (vButton == mChangeParametersButton) {
                PFTracking3D.this.getUserDefinedParams();
            }

            //
            // Calculate Button
            //
            if (vButton == mCalcFromHereButton) {
                PFTracking3D.this.calcFromHereButtonPressed();
            }

            //
            // Save init button
            //
            if (vButton == mSaveInitButton) {
                PFTracking3D.this.saveInitButtonPressed();
            }

            //
            // Save results button
            //
            if (vButton == mSaveCorrectionButton) {
                PFTracking3D.this.saveCorrectionButtonPressed();
            }

            //
            // Clear Button
            //
            if (vButton == mDeleteAllButton) {
                PFTracking3D.this.deleteAllButtonPressed();
            }

            //
            // initialze with mouse button
            //
            if (vButton == mMouseInitializationButton) {
                PFTracking3D.this.initializeWithMouseButtonPressed();
            }

            // generate an updated view with the ImagePlus in this window according to the new filter
//			generateView(this.imp);
        }

        /**
         * Defines the action taken upon an
         * <code>MouseEvent</code> triggered by left-clicking the mouse anywhere
         * in this
         * <code>TrajectoryStackWindow</code>
         *
         * @see
         * java.awt.event.MouseListener#mousePressed(java.awt.event.MouseEvent)
         */
        public synchronized void mousePressed(MouseEvent aE) {
            PFTracking3D.this.mousePressed(this.ic.offScreenX(aE.getPoint().x), this.ic.offScreenY(aE.getPoint().y));
        }

        public void mouseReleased(MouseEvent aE) {
            PFTracking3D.this.mouseReleased(this.ic.offScreenX(aE.getPoint().x), this.ic.offScreenY(aE.getPoint().y));
        }

        public void mouseClicked(MouseEvent aE) {
            PFTracking3D.this.mouseClicked(this.ic.offScreenX(aE.getPoint().x), this.ic.offScreenY(aE.getPoint().y));
        }

        public void mouseEntered(MouseEvent aE) {
            PFTracking3D.this.mouseEntered(aE);
        }

        public void mouseExited(MouseEvent aE) {
            PFTracking3D.this.mouseExited(aE);
        }
    } // CustomStackWindow inner class

    public STATE_OF_FILTER getMStateOfFilter() {
        return mStateOfFilter;
    }

    public void setMStateOfFilter(STATE_OF_FILTER stateOfFilter) {
        mStateOfFilter = stateOfFilter;
    }

    public int getMNbThreads() {
        return mNbThreads;
    }

    public void setMNbThreads(int nbThreads) {
        mNbThreads = nbThreads;
    }

    public int getMNbParticles() {
        return mNbParticles;
    }

    public void setMNbParticles(int nbParticles) {
        mNbParticles = nbParticles;
    }

    public int getMRepSteps() {
        return mRepSteps;
    }

    public void setMRepSteps(int repSteps) {
        mRepSteps = repSteps;
    }

    public int getMInitParticleFilterIterations() {
        return mInitRWIterations;
    }

    public void setMInitParticleFilterIterations(int initParticleFilterIterations) {
        mInitRWIterations = initParticleFilterIterations;
    }

    public int getMResamplingThreshold() {
        return mResamplingThreshold;
    }

    public void setMResamplingThreshold(int resamplingThreshold) {
        mResamplingThreshold = resamplingThreshold;
    }

    public float getMBackground() {
        return mBackground;
    }

    public void setMBackground(float background) {
        mBackground = background;
    }

    abstract public float[] getMSigmaOfRandomWalk();

    abstract public String[] getMDimensionsDescription();

    abstract public boolean getMDoPrecisionOptimization();

    public void setMDimensionsDescription(String[] dimensionsDescription) {
        mDimensionsDescription = dimensionsDescription;
    }

    public long getMSeed() {
        return mSeed;
    }

    public void setMSeed(long seed) {
        mSeed = seed;
    }

    public boolean isMDoResampling() {
        return mDoResampling;
    }

    public void setMDoResampling(boolean doResampling) {
        mDoResampling = doResampling;
    }

    public boolean isMDoPrintStates() {
        return mDoPrintStates;
    }

    public void setMDoPrintStates(boolean doPrintStates) {
        mDoPrintStates = doPrintStates;
    }

    public boolean isMDoMonitorIdealImage() {
        return mDoMonitorIdealImage;
    }

    public void setMDoMonitorIdealImage(boolean doMonitorIdealImage) {
        mDoMonitorIdealImage = doMonitorIdealImage;
    }

    public boolean isMDoMonitorParticles() {
        return mDoMonitorParticles;
    }

    public void setMDoMonitorParticles(boolean doMonitorParticles) {
        mDoMonitorParticles = doMonitorParticles;
    }

    public ImageStack getMIdealImageMonitorProcessor() {
        return mIdealImageMonitorStack;
    }

    public void setMIdealImageMonitorProcessor(ImageStack idealImageMonitorProcessor) {
        mIdealImageMonitorStack = idealImageMonitorProcessor;
    }

    public Vector<Vector<Vector<float[]>>> getMParticleMonitor() {
        return mParticleMonitor;
    }

    public void setMParticleMonitor(Vector<Vector<Vector<float[]>>> particleMonitor) {
        mParticleMonitor = particleMonitor;
    }

    public int getMFrameOfInitialization() {
        return mFrameOfInitialization;
    }

    public void setMFrameOfInitialization(int frameOfInitialization) {
        mFrameOfInitialization = frameOfInitialization;
    }

    public int getMWidth() {
        return mWidth;
    }

    public void setMWidth(int width) {
        mWidth = width;
    }

    public ImagePlus getMZProjectedImagePlus() {
        return mZProjectedImagePlus;
    }

    public void setMZProjectedImagePlus(ImagePlus projectedImagePlus) {
        mZProjectedImagePlus = projectedImagePlus;
    }

    public Vector<Vector<float[]>> getMParticles() {
        return mParticles;
    }

    public void setMParticles(Vector<Vector<float[]>> particles) {
        mParticles = particles;
    }

    public static String getRESULT_FILE_SUFFIX() {
        return RESULT_FILE_SUFFIX;
    }

    public static String getINIT_FILE_SUFFIX() {
        return INIT_FILE_SUFFIX;
    }

    public int getMHeight() {
        return mHeight;
    }

    public int getMNSlices() {
        return mNSlices;
    }

    public int getMNFrames() {
        return mNFrames;
    }

    public ImagePlus getMOriginalImagePlus() {
        return mOriginalImagePlus;
    }

    public Vector<Vector<float[]>> getMStateVectorsMemory() {
        return mStateVectorsMemory;
    }

    public float getMSigmaPSFxy() {
        return mSigmaPSFxy;
    }

    public void setMSigmaPSFxy(float sigmaPSFxy) {
        mSigmaPSFxy = sigmaPSFxy;
    }

    public float getMSigmaPSFz() {
        return mSigmaPSFz;
    }

    public void setMSigmaPSFz(float sigmaPSFz) {
        mSigmaPSFz = sigmaPSFz;
    }

    public float getMNA() {
        return mNA;
    }

    public void setMNA(float mna) {
        mNA = mna;
    }

    public float getMn() {
        return mn;
    }

    public void setMn(float mn) {
        this.mn = mn;
    }

    public float getPixelWidthInNm() {
        return (float) mOriginalImagePlus.getCalibration().pixelWidth;
    }

    public float getPixelDepthInNm() {
        return (float) mOriginalImagePlus.getCalibration().pixelDepth;
    }

    public class Line3D {

        Point3D mA, mB;

        public Line3D() {
            this(new Point3D(), new Point3D());
        }

        public Line3D(Point3D aStartPoint, Point3D aEndPoint) {
            mA = aStartPoint;
            mB = aEndPoint;
        }

        public Point3D getMA() {
            return mA;
        }

        public void setMA(Point3D ma) {
            mA = ma;
        }

        public Point3D getMB() {
            return mB;
        }

        public void setMB(Point3D mb) {
            mB = mb;
        }

        public float getDistanceToLine(Point3D aP) {
            //|(p-a)x(b-a)| / |b-a| , u = b-a
            Point3D vU = mB.clone().subtract(mA);
            return aP.clone().subtract(mA).cross(vU).getLength() / vU.getLength();
        }

        public float getDistanceToSegment(Point3D aPoint) {
            if (aPoint.clone().subtract(mA).normalize().scalarProduct(mB.clone().subtract(mA).normalize()) < 0) {
                return aPoint.clone().subtract(mA).getLength();
            }
            if (aPoint.clone().subtract(mB).normalize().scalarProduct(mA.clone().subtract(mB).normalize()) > 0) {
                return aPoint.clone().subtract(mB).getLength();
            }
            return getDistanceToLine(aPoint);
        }

        /**
         *
         * @return a Line from the closest point to (0,0,0) to the point most
         * far away.
         */
        public Line3D getBoundingBox() {
            return new Line3D(
                    new Point3D(Math.min(mA.mX, mB.mX), Math.min(mA.mY, mB.mY), Math.min(mA.mZ, mB.mZ)),
                    new Point3D(Math.max(mA.mX, mB.mX), Math.max(mA.mY, mB.mY), Math.max(mA.mZ, mB.mZ)));
        }
    }

    public class Point3D {

        public float mX, mY, mZ;

        public Point3D() {
            this(0f, 0f, 0f);
        }

        public Point3D(float aX, float aY, float aZ) {
            mX = aX;
            mY = aY;
            mZ = aZ;
        }

        public Point3D clone() {
            return new Point3D(mX, mY, mZ);
        }

        public Point3D add(Point3D aB) {
            mX += aB.mX;
            mY += aB.mY;
            mZ += aB.mZ;
            return this;
        }

        /**
         * this = this - aB
         *
         * @param aB
         */
        public Point3D subtract(Point3D aB) {
            mX -= aB.mX;
            mY -= aB.mY;
            mZ -= aB.mZ;
            return this;
        }

        public float scalarProduct(Point3D aB) {
            return mX * aB.mX + mY * aB.mY + mZ * aB.mZ;
        }

        public Point3D cross(Point3D aB) {
            mX = mY * aB.mZ - mZ * aB.mY;
            mY = mZ * aB.mX - mX * aB.mZ;
            mZ = mX * aB.mY - mY * aB.mX;
            return this;
        }

        public float getLength() {
            return (float) Math.sqrt(scalarProduct(this));
        }

        public Point3D normalize() {
            float vL = getLength();
            mX /= vL;
            mY /= vL;
            mZ /= vL;
            return this;
        }

        public float getMX() {
            return mX;
        }

        public void setMX(float mx) {
            mX = mx;
        }

        public float getMY() {
            return mY;
        }

        public void setMY(float my) {
            mY = my;
        }

        public float getMZ() {
            return mZ;
        }

        public void setMZ(float mz) {
            mZ = mz;
        }
    }
}
