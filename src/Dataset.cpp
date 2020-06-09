#include "def.h"

class Dataset
{
    // contigency table index translation
    // note that 3 or 0b11 is not a valid genotype and should not be considered in Gini Computation.
    uint64 CaseIndex(varIdx v, sampleIdx s) { return ((v * numByteCase) + s); } // get the byte index of sample in data
    uint64 CtrlIndex(varIdx v, sampleIdx s) { return ((v * numByteCtrl) + s); } // get the byte index of sample in data

public:
    uint32 order;

    bool *labels;

    sampleIdx numSamples;

    sampleIdx numCase;
    sampleIdx numCtrl;

    uint32 numWordCase; // number of machine word used to store Case data (each sample is a byte)
    uint32 numWordCtrl; // number of machine word used to store Ctrl data (each sample is a byte)

    uint32 numByteCase; // numWordCase * sizeof(word)
    uint32 numByteCtrl; // numWordCtrl * sizeof(word)

    uint8 *byteCase[MAX_ORDER]; // byte pointer to store genotype data and their shifted version
    uint8 *byteCtrl[MAX_ORDER]; // byte pointer to store genotype data and their shifted version

    word *wordCase[MAX_ORDER]; // the wrod pointer to byteCase
    word *wordCtrl[MAX_ORDER]; // the wrod pointer to byteCtrl

    varIdx numVariables;
    char **nameVariable;

    double setBeta; // beta of the original set

    sampleIdx *contingency_table; // should be small enough to remain in cache

    void FreeMemory()
    {
        delete[] labels;

        for (uint32 i = 0; i < numVariables; i++)
            delete[] nameVariable[i];
        delete nameVariable;

        for (uint32 i = 0; i < order; i++)
        {
            delete[] wordCase[i];
            delete[] wordCtrl[i];
        }
    }

    void prepareDataset(std::vector<bool> &sampleClasses, std::vector<std::string> &snpLabels, std::vector<std::vector<uint8_t>> &genotypes)
    {
        numSamples = sampleClasses.size();
        if (numSamples >= pow(2, sizeof(sampleIdx) * 8))
            ERROR("Change sampleIdx type to support the number of samples exist in dataset");

        numVariables = snpLabels.size();
        numCase = std::count(sampleClasses.begin(), sampleClasses.end(), true);
        numCtrl = numSamples - numCase;

        //prepare sample class labels array
        labels = new bool[numSamples];
        NULL_CHECK(labels);
        std::copy(sampleClasses.begin(), sampleClasses.end(), labels);

        // find number of word and byte per variable in Case and Ctrl
        numWordCase = numCase / byte_in_word;
        numWordCtrl = numCtrl / byte_in_word;

        if (numCase % byte_in_word)
            numWordCase++;
        if (numCtrl % byte_in_word)
            numWordCtrl++;

        numByteCase = numWordCase * sizeof(word);
        numByteCtrl = numWordCtrl * sizeof(word);

        // allocate memory
        wordCase[0] = new word[numVariables * numWordCase];
        wordCtrl[0] = new word[numVariables * numWordCtrl];

        NULL_CHECK(wordCase[0]);
        NULL_CHECK(wordCtrl[0]);

        // convert to byte address
        byteCase[0] = (uint8_t *)wordCase[0];
        byteCtrl[0] = (uint8_t *)wordCtrl[0];

        //prepare the array of char*s for the names of each SNP
        nameVariable = new char *[numVariables];
        NULL_CHECK(nameVariable);
        for (int i = 0; i < snpLabels.size(); ++i)
        {
            nameVariable[i] = new char[snpLabels[i].size() + 1];
            strcpy(nameVariable[i], snpLabels[i].c_str());
        }

        //prepare the arrays of the variants of each sample
        for (int SNP = 0; SNP < genotypes.size(); SNP++)
        {
            uint32_t idxCase = 0;
            uint32_t idxCtrl = 0;

            for (sampleIdx i = 0; i < numSamples; i++)
            {
                if (sampleClasses[i])
                {
                    byteCase[0][CaseIndex(SNP, idxCase)] = genotypes[SNP][i];
                    idxCase++;
                }
                else
                {
                    byteCtrl[0][CtrlIndex(SNP, idxCtrl)] = genotypes[SNP][i];
                    idxCtrl++;
                }
            }
        }

        std::cout << "There are " << snpLabels.size() << " SNPs" << std::endl;
        std::cout << "There are " << sampleClasses.size() << " samples" << std::endl;
        std::cout << "There are " << numCase << " Cases" << std::endl;
        std::cout << "There are " << numCtrl << " Controls" << std::endl;
    }

    void ReadDatasetBfile(const char *fn)
    {
        std::string filename = std::string(fn);
        //if the file ends with .bed strip it
        if (filename.size() > 4 && filename.compare(filename.size() - 4, 4, ".bed") == 0)
            filename.erase(filename.end() - 4, filename.end());

        std::cout << "loading dataset " << filename << ".bed, "
                  << filename << ".bim, " << filename << ".fam" << std::endl;

        std::vector<bool> sampleClasses;
        std::vector<std::string> snpLabels;
        std::vector<std::vector<uint8_t>> genotypes;
        //get variant labels
        std::ifstream variantFile(filename + ".bim");
        if (!variantFile.is_open())
            ERROR("could not open associated .bim file");
        std::string line;
        while (std::getline(variantFile, line))
        {
            //find delimeter
            char delim = ' ';
            if (std::count(line.begin(), line.end(), '\t') == 5)
                delim = '\t';
            line = line.substr(line.find(delim) + 1); //remove first field
            snpLabels.push_back(line.substr(0, line.find(delim)));
        }
        variantFile.close();

        //get sample classes
        std::ifstream sampleFile(filename + ".fam");
        if (!sampleFile.is_open())
            ERROR("could not open associated .fam file");
        while (std::getline(sampleFile, line))
        {
            std::string classLabel = line.substr(line.find_last_not_of(" \v\f\t\r\n"), 1);
            if (classLabel == "1")
                sampleClasses.push_back(false);
            else if (classLabel == "2")
                sampleClasses.push_back(true);
            else
                ERROR("Phenotype type value not control (1) or case (2) in .fam file")
        }
        sampleFile.close();

        //get sample variants
        std::ifstream bedFile(filename + ".bed", std::ios::binary);
        if (!bedFile.is_open())
            ERROR("could not open associated .bed file");
        uint32_t chunkSize = sampleClasses.size() / 4;
        if (sampleClasses.size() % 4 != 0)
            chunkSize++;
        char *buffer = new char[chunkSize];
        genotypes.reserve(snpLabels.size());
        //read and discard 3 magic bytes at the start of the file
        bedFile.read(buffer, 3);
        for (int i = 0; i < snpLabels.size(); ++i)
        {
            bedFile.read(buffer, chunkSize);
            if (!bedFile.good())
                ERROR("bed file does not contain enough bytes")
            genotypes.push_back(std::vector<uint8_t>());
            std::vector<uint8_t> &gs = genotypes.back();
            gs.reserve(chunkSize * 4);
            for (int j = 0; j < chunkSize; ++j)
            {
                //genotypes stored as two bit value
                //00 homozygous minor allele recorded as 2
                //01 missing value recorded as homozygous major allele recorded as 0
                //10 heterogenous recorded as 1
                //11 homozygous major allele recorded as 0
                static const uint8_t convert[4] = {2, 0, 1, 0};
                gs.push_back(convert[buffer[j] & 3]);
                gs.push_back(convert[(buffer[j] >> 2) & 3]);
                gs.push_back(convert[(buffer[j] >> 4) & 3]);
                gs.push_back(convert[(buffer[j] >> 6) & 3]);
            }
            //discard the last samples that were created to pad out the 8 bits
            gs.resize(sampleClasses.size());
        }

        prepareDataset(sampleClasses, snpLabels, genotypes);
        return;
    }

    void ComputeSetBeta()
    {
        numCase = 0;
        numCtrl = 0;
        for (sampleIdx i = 0; i < numSamples; i++)
            if (labels[i])
                numCase++;
            else
                numCtrl++;
        setBeta = P2((double)numCase / numSamples) + P2((double)numCtrl / numSamples);
        printf("\n Purity of the whole dataset (B_0) is %f (baseline for Beta)", setBeta);
    }

    void Shift()
    {
        for (uint32 d = 1; d < order; d++)
        {
            // allocate memory
            wordCase[d] = new word[numVariables * numWordCase];
            wordCtrl[d] = new word[numVariables * numWordCtrl];
            NULL_CHECK(wordCase[d]);
            NULL_CHECK(wordCtrl[d]);
            // convert to byte address
            byteCase[d] = (uint8 *)wordCase[d];
            byteCtrl[d] = (uint8 *)wordCtrl[d];

            // shoft and copy
            for (uint32 i = 0; i < (numVariables * numWordCase); i++)
                wordCase[d][i] = wordCase[d - 1][i] << 2;
            for (uint32 i = 0; i < (numVariables * numWordCtrl); i++)
                wordCtrl[d][i] = wordCtrl[d - 1][i] << 2;

            printf("\n Shift dataset by %u bits compeleted", d * 2);
        }
    }

    void Init(ARGS args)
    {
        order = args.maxOrder;
        ComputeSetBeta();
        Shift();
    }

    word *GetVarCase(uint32 o, varIdx vi)
    {
        return &wordCase[o][vi * numWordCase];
    }

    word *GetVarCtrl(uint32 o, varIdx vi)
    {
        return &wordCtrl[o][vi * numWordCtrl];
    }
};
