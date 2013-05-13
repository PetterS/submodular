/*//////////////////////////////////////////////////////////////////////////////////////////////////
///  Generic Cuts Software
///  Version 1.0
////////////////////////////////////////////////////////////////////////////////////////////////////

Copyright 2013 Chetan Arora.
This software can be used for research purposes only. This software or its derivatives must not be publicly distributed without a prior consent from the author (Chetan Arora).

THIS SOFTWARE IS PROVIDED "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

For the latest version, check: http://www.cs.huji.ac.il/~chetan

The software implements the technique described in the following paper:

Chetan Arora, Subhashis Banerjee, Prem Kalra and S.N. Maheshwari. "Generic Cuts: An efficient optimal algorithm for submodular MRF-MAP problems with higher order cliques". In ECCV 2012.

////////////////////////////////////////////////////////////////////////////////////////////////////

This software is implemented so that it can be used most conveniently in combination with the OpenCV Software. 

The code was tested on Microsoft Visual C++ 2008, and 2010 Professional with OpenCV 2.4.3 on Windows 7, 64 Bit. 

////////////////////////////////////////////////////////////////////////////////////////////////////

The software implements the Augmenting Path style algorithm as described in the ECCV paper as well as the Push Relabel style for finding max flow in the gadget based flow graph. Abstract class GadgetGraph in the util folder gives a common interface to the two implemented styles as well as future planned releases.

The mainfile.cpp included in the software describes a sample way to to use the software. 