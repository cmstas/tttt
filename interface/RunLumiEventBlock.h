#ifndef RUNEVENTLUMIBLOCK_H
#define RUNEVENTLUMIBLOCK_H


typedef unsigned int RunNumber_t;
typedef unsigned int LuminosityBlock_t;
typedef unsigned long long EventNumber_t;

#define RUN_VARIABLE \
RUNLUMIEVENT_VARIABLE(RunNumber_t, RunNumber, run)
#define LUMI_VARIABLE \
RUNLUMIEVENT_VARIABLE(LuminosityBlock_t, LuminosityBlock, luminosityBlock)
#define EVENT_VARIABLE \
RUNLUMIEVENT_VARIABLE(EventNumber_t, EventNumber, event)

#define RUNLUMI_VARIABLES \
RUN_VARIABLE \
LUMI_VARIABLE

#define RUNLUMIEVENT_VARIABLES \
RUNLUMI_VARIABLES \
EVENT_VARIABLE


class RunLumiEventBlock{
protected:
#define RUNLUMIEVENT_VARIABLE(TYPE, NAME, NANONAME) TYPE NAME;
  RUNLUMIEVENT_VARIABLES;
#undef RUNLUMIEVENT_VARIABLE

public:
  RunLumiEventBlock();
  RunLumiEventBlock(unsigned int RunNumber_, unsigned int LuminosityBlock_, unsigned long long EventNumber_);
  RunLumiEventBlock(RunLumiEventBlock const& other);

  bool operator==(RunLumiEventBlock const& other) const;
};


#endif
