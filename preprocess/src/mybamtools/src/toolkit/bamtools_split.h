// ***************************************************************************
// bamtools_split.h (c) 2010 Derek Barnett, Erik Garrison
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ---------------------------------------------------------------------------
// Last modified: 18 September 2010 (DB)
// ---------------------------------------------------------------------------
// 
// ***************************************************************************

#ifndef BAMTOOLS_SPLIT_H
#define BAMTOOLS_SPLIT_H

#include "bamtools_tool.h"

namespace BamTools {
  
class SplitTool : public AbstractTool {
  
    public:
        SplitTool(void);
        ~SplitTool(void);
  
    public:
        int Help(void);
        int Run(int argc, char* argv[]); 
        
    private:
        struct SplitSettings;
        SplitSettings* m_settings;
        
        struct SplitToolPrivate;
        SplitToolPrivate* m_impl;
};
  
} // namespace BamTools

#endif // BAMTOOLS_SPLIT_H