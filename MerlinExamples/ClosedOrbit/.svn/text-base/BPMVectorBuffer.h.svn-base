#include "AcceleratorModel/ActiveMonitors/BPM.h"

class BPMVectorBuffer : public BPM::Buffer
{
public:
	BPMVectorBuffer () {};

	vector<BPM::Data> BPMReading;

	void Record (const BPM& aBPM, const BPM::Data& data) {
		BPMReading.push_back(data);
	};

	void Clear () {
		BPMReading.clear();
	};
};

