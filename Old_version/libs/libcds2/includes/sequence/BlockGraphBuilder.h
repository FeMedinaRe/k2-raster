//
// Created by alberto on 4/8/15.
//
#include <sequence/SequenceBuilder.h>
#include <sequence/BlockGraphUtils.h>
#include <sequence/BlockGraph.h>


#ifndef _ALL_SEQUENCEBUILDERSEQUENCEBLOCKGRAPH_H_
#define _ALL_SEQUENCEBUILDERSEQUENCEBLOCKGRAPH_H_

namespace cds_static {
    class BlockGraphBuilder {
    public:
        BlockGraphBuilder(uint nLevels, uint blockLengthAtLeaves) {
            this->nLevels = nLevels;
            this->blockLengthAtLeaves = blockLengthAtLeaves;
            gs=2;
            maxV=DEFAULT_MAX_SIGMA;
            inputBufferSize=DEFAULT_BUFFER_READ_SIZE;
            outputBufferSize=DEFAULT_BUFFER_WRITE_SIZE;
            bsb=nullptr;
        }

        virtual BlockGraphBuilder& setMarksBuilder(BitSequenceBuilder *_bsb){
            this->bsb = _bsb;
            return *this;
        }
        /**
        * @param outStream: stream to write out the data structure. It must be opened if it is
        not empty. The data structure is written level by level, freeing the space the data
        structures take after written them.
        */
        virtual BlockGraphBuilder& setOutStream(string &_outFile){
            this->outFile = _outFile;
            return *this;
        }
        virtual BlockGraphBuilder& setMaxV(uint _maxv){
            this->maxV = _maxv;
            return *this;
        }
        virtual BlockGraphBuilder& setGroupSize(uint _gs){
            this->gs = _gs;
            return *this;
        }
        virtual BlockGraphBuilder& setBufferInputReadSize(ulong _inputBufferSize){
            this->inputBufferSize = _inputBufferSize;
            return *this;
        }
        virtual BlockGraphBuilder& setBufferOutputReadSize(ulong _outputBufferSize){
            this->outputBufferSize = _outputBufferSize;
            return *this;
        }
        virtual BlockGraphBuilder& setBufferInName(string && _tmpFileNameBufferIn){
            this->tmpFileNameBufferIn = std::move(_tmpFileNameBufferIn);
            return *this;
        }
        virtual BlockGraphBuilder& setBufferOutName(string &&_tmpFileNameBufferOut){
            this->tmpFileNameBufferOut = std::move(_tmpFileNameBufferOut);
            return *this;
        }

        virtual ~BlockGraphBuilder() {
        }

        BlockGraph<uint,Sequence> *build(uint *seq, size_t len) {
            return nullptr;
        }

        BlockGraph<uint,Sequence> *build(string &file) {
            StrategySplitBlock *st = new StrategySplitBlockBottomUp(nLevels, blockLengthAtLeaves);
            //create the block graph (does not run the construction algorithm)
            return new BlockGraph<uint, Sequence>(file,st,bsb,outFile,maxV,inputBufferSize,
                    outputBufferSize,tmpFileNameBufferIn,tmpFileNameBufferOut,gs);
        }

    protected:
        uint gs;
        string outFile;
        uint maxV;
        ulong inputBufferSize;
        ulong outputBufferSize;
        string tmpFileNameBufferIn;
        string tmpFileNameBufferOut;
        uint nLevels;
        uint blockLengthAtLeaves;
        BitSequenceBuilder *bsb;
    };

    class SequenceBlockGraphBuilder : public BlockGraphBuilder,public SequenceBuilder{
    public:

        SequenceBlockGraphBuilder(uint nLevels, uint blockLengthAtLeaves):BlockGraphBuilder(nLevels, blockLengthAtLeaves) {
            seq_builder = nullptr;
        }
        SequenceBlockGraphBuilder& setSequenceBuilder(SequenceBuilder *sb){
            this->seq_builder = sb;
            return *this;
        }

        virtual ~SequenceBlockGraphBuilder(){}

        virtual Sequence *build(uint *seq, size_t len) {
            TemporaryFile tmpFile;
            ofstream fp(tmpFile.getFileName());
            fp.write((char*)seq,sizeof(uint)*len);
            fp.close();
            string file_input(tmpFile.getFileName());
            StrategySplitBlock *st = new StrategySplitBlockBottomUp(nLevels, blockLengthAtLeaves);
            //create the block graph (does not run the construction algorithm)
            return new SequenceBlockGraph<uint, Sequence>(file_input,st,bsb,outFile,maxV,inputBufferSize,
                                                          outputBufferSize,tmpFileNameBufferIn,tmpFileNameBufferOut,gs, seq_builder);
        }

        virtual Sequence *build(const Array &tmp) {
            return nullptr;
        }
        virtual Sequence *build(string &file){
            StrategySplitBlock *st = new StrategySplitBlockBottomUp(nLevels, blockLengthAtLeaves);
            //create the block graph (does not run the construction algorithm)
            return new SequenceBlockGraph<uint, Sequence>(file,st,bsb,outFile,maxV,inputBufferSize,
                    outputBufferSize,tmpFileNameBufferIn,tmpFileNameBufferOut,gs, seq_builder);

        }

    protected:

        SequenceBuilder *seq_builder;
    };
};
#endif //_ALL_SEQUENCEBUILDERSEQUENCEBLOCKGRAPH_H_
