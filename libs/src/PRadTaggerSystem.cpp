//============================================================================//
// Tagger system class                                                        //
// Decode tagger information from raw data                                    //
//                                                                            //
// Chao Peng                                                                  //
// 11/20/2016                                                                 //
//============================================================================//

#include "PRadTaggerSystem.h"
#include "TH2I.h"



//============================================================================//
// Constructors, Destructor, Assignment Operators                             //
//============================================================================//

// constructor
PRadTaggerSystem::PRadTaggerSystem()
{
    hist_E = new TH2I("Tagger E", "Tagger E counter", 2000, 0, 20000, 384, 0, 383);
    hist_T = new TH2I("Tagger T", "Tagger T counter", 2000, 0, 20000, 128, 0, 127);
}

// copy/move constructors
PRadTaggerSystem::PRadTaggerSystem(const PRadTaggerSystem &that)
{
    hist_E = new TH2I(*that.hist_E);
    hist_T = new TH2I(*that.hist_T);
}

PRadTaggerSystem::PRadTaggerSystem(PRadTaggerSystem &&that)
{
    hist_E = that.hist_E;
    that.hist_E = nullptr;
    hist_T = that.hist_T;
    that.hist_T = nullptr;
}

// destructor
PRadTaggerSystem::~PRadTaggerSystem()
{
    delete hist_E;
    delete hist_T;
}

// copy/move assignment operators
PRadTaggerSystem &PRadTaggerSystem::operator =(const PRadTaggerSystem &rhs)
{
    if(this == &rhs)
        return *this;

    PRadTaggerSystem that(rhs);
    *this = std::move(that);
    return *this;
}

PRadTaggerSystem &PRadTaggerSystem::operator =(PRadTaggerSystem &&rhs)
{
    if(this == &rhs)
        return *this;

    delete hist_E;
    delete hist_T;

    hist_E = rhs.hist_E;
    hist_T = rhs.hist_T;

    return *this;
}



//============================================================================//
// Public Member Functions                                                    //
//============================================================================//

// reset current data and hists
void PRadTaggerSystem::Reset()
{
    hist_E->Reset();
    hist_T->Reset();
}

// feed tagger hits to event data
void PRadTaggerSystem::FeedTaggerHits(const TDCV1190Data &tdcData, EventData &event)
{
    // E channel
    if(tdcData.addr.slot == 3 || tdcData.addr.slot == 5 || tdcData.addr.slot == 7)
    {
        int e_ch = tdcData.addr.channel + (tdcData.addr.slot - 3)*64 + TAGGER_CHANID;
        // E Channel 30000 + channel
        event.add_tdc(TDC_Data(e_ch, tdcData.val));
    }
    // T Channel
    if(tdcData.addr.slot == 14)
    {
        int t_lr = tdcData.addr.channel/64;
        int t_ch = tdcData.addr.channel%64;
        if(t_ch > 31)
            t_ch = 32 + (t_ch + 16)%32;
        else
            t_ch = (t_ch + 16)%32;
        t_ch += t_lr*64;
        event.add_tdc(TDC_Data(t_ch + TAGGER_CHANID + TAGGER_T_CHANID, tdcData.val));
    }
}

// fill tagger hits to histograms
void PRadTaggerSystem::FillHists(const EventData &event)
{
    for(auto &tdc : event.get_tdc_data())
    {
        if(tdc.channel_id < TAGGER_CHANID)
            continue;
        int id = tdc.channel_id - TAGGER_CHANID;
        if(id >= TAGGER_T_CHANID)
            hist_T->Fill(tdc.value, id - TAGGER_T_CHANID);
        else
            hist_E->Fill(tdc.value, id);
    }
}
