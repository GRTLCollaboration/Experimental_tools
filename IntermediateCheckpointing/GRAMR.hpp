#ifndef GRAMR_HPP_
#define GRAMR_HPP_

#include "AMR.H"

class GRAMR : public AMR
{
  public:
    /// Set up for a restart from an intermediate checkpoint file.
    /// Returns true if we can continue
    bool prepareIntermediateRestart();

  protected:
    /// This function advances levels with different times until they are at the
    /// same time (e.g. after a restart from an intermediate checkpoint file)
    void bringLevelsIntoSync(std::vector<int> &a_steps_until_sync,
                             std::vector<bool> &a_step_finished);

    void writeIntermediateCheckpointFile(
        const std::string additional_prefix = "INTERMEDIATE_");

    bool isSteadyState();
};

#endif
