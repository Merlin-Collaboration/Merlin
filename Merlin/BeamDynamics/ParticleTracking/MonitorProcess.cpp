#include "MonitorProcess.h"
#include "AcceleratorModel/AcceleratorComponent.h"
#include <fstream>

MonitorProcess::MonitorProcess(const string& aID ,  int prio, const string& prefix ):BunchProcess(aID,prio)
{
        cout << "MonitorProcess()" << endl;
        active = false;
        file_prefix = prefix;
}

void MonitorProcess::SetPrefix(const string& prefix)
{
        file_prefix = prefix;
}

void MonitorProcess::AddElement(const string e){
        dump_at_elements.push_back(e);
}


void MonitorProcess::InitialiseProcess (Bunch& this_bunch){
        bunch = &this_bunch;
//      cout << "MonitorProcess::InitialiseProcess()" << endl;
}
void MonitorProcess::DoProcess (const double ds){
        string filename;
        filename = file_prefix + currentComponent->GetName();

        cout << "MonitorProcess::DoProcess(): " << filename << endl;
        ofstream out_file(filename.c_str());
        bunch->Output(out_file);
        out_file.close();
}

double MonitorProcess::GetMaxAllowedStepSize() const{
        return 1000;
}

void MonitorProcess::SetCurrentComponent (AcceleratorComponent& component)
{
        currentComponent = &component;
        vector<string>::iterator result;
        result = find(dump_at_elements.begin(), dump_at_elements.end(), component.GetName());

        if( result == dump_at_elements.end() )
                active = false;
        else
                active = true;
}
