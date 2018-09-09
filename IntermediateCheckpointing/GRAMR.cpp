#include "GRAMR.hpp"
#include "GRAMRLevel.hpp"

bool GRAMR::prepareIntermediateRestart()
{
    std::vector<int> steps_until_sync(m_finest_level + 1);
    std::vector<bool> step_finished(m_finest_level + 1);

    steps_until_sync[0] = 0;
    for (int ilev = 1; ilev <= m_finest_level; ++ilev)
    {
        Real time = m_amrlevels[ilev]->time();
        Real coarser_time = m_amrlevels[ilev - 1]->time();
        Real dt = m_amrlevels[ilev]->dt();

        // dt*0.1 prevents problems with floating point arithmetic
        steps_until_sync[ilev] = (coarser_time - time + dt * 0.1) / dt;
    }

    // Prepare the step_finished vector and output
    pout() << "level   time   steps_left   step_finished\n"; // Write table
    for (int ilev = 0; ilev <= m_finest_level; ++ilev)
    {
        GRAMRLevel *gr_level_ptr = gr_cast(m_amrlevels[ilev]);
        step_finished[ilev] = gr_level_ptr->step_finished();
        pout() << "  " << ilev << "    " << m_amrlevels[ilev]->time()
               << "       " << steps_until_sync[ilev] << "              "
               << step_finished[ilev] << endl;
    }

    bool can_continue = true;
    // Only do stuff if the coarsest level hasn't finished the step
    if (step_finished[0] == false)
    {
        pout() << "Synchronising levels..." << endl;
        bringLevelsIntoSync(steps_until_sync, step_finished);
        // Update the counters
        ++m_cur_step;
        m_cur_time += m_dt_base;

        // If steady state has been reached we shouldn't continue
        can_continue = !isSteadyState();
        // Check whether the levels are now synchronised.
        if (steps_until_sync != std::vector<int>(steps_until_sync.size(), 0))
        {
            MayDay::Error(
                "Levels failed to synchronise. Level times in checkpoint "
                "might be incompatible with subcycling procedure.");
            can_continue = false;
        }
    }
    else
    {
        // If the coarsest level is finished all other level must be in sync
        CH_assert(steps_until_sync == std::vector<int>(m_finest_level + 1, 0));
    }
    return can_continue;
}

void GRAMR::bringLevelsIntoSync(std::vector<int> &a_steps_until_sync,
                                std::vector<bool> &a_step_finished)
{
    CH_assert(isDefined());
    CH_assert(isSetUp());
    CH_assert(a_steps_until_sync.size() == m_finest_level + 1);
    // If the coarsest step is done you shouldn't be calling this function
    CH_assert(!a_step_finished[0]);
    // If the finest level is in the middle of a step something is wrong...
    CH_assert(a_step_finished[m_finest_level]);
    // The corsest level shouldn't have to sync with anything
    CH_assert(a_steps_until_sync[0] == 0);

    // Traverse levels from fine to coarse and synchronise the levels
    // Coarser levels will subcylce finer levels if necessary
    for (int ilev = m_finest_level; ilev > 0; --ilev)
    {
        bool timeBoundary = true;
        int steps_left = a_steps_until_sync[ilev];
        if (steps_left > 0)
        {
            // If we still need to do steps the level above can't be finished
            CH_assert(a_step_finished[ilev - 1] != true);

            a_steps_until_sync[ilev] = timeStep(ilev, steps_left, timeBoundary);
            m_amrlevels[ilev - 1]->postTimeStep();
        }
        else
        {
            // If the level didn't have to be advanced but step of the coarser
            // level isn't finished yet the only thing left to do is to call
            // postTimeStep
            if (!a_step_finished[ilev - 1])
            {
                m_amrlevels[ilev - 1]->postTimeStep();
            }
        }
    }
    assignDt(); // In AMR.cpp this is called after the coarsest level step
}

void GRAMR::writeIntermediateCheckpointFile(const std::string additional_prefix)
{
    std::string tmp = m_checkpointfile_prefix;
    m_checkpointfile_prefix += additional_prefix;
    writeCheckpointFile();
    m_checkpointfile_prefix = tmp;
}

bool GRAMR::isSteadyState()
{
    // TODO: This code is copied verbatim from AMR.cpp ... ideally amr.run in
    // Chombo should be broken up.
    bool steadyStateStop = false;
    if (m_checkForSteadyState)
    {
        steadyStateStop = m_amrlevels[0]->convergedToSteadyState();
        for (int level = 1; level <= m_finest_level; ++level)
        {
            steadyStateStop = ((steadyStateStop) &&
                               (m_amrlevels[level]->convergedToSteadyState()));
        }
    }
    if (steadyStateStop)
    {
        pout() << "GRAMR: steady state was reached.\n";
        if (m_plot_interval > 0)
        {
            pout() << "Writing plot file...\n";
            writePlotFile();
        }
        if ((m_checkpoint_interval > 0) && (m_lastcheck_step != m_cur_step))
        {
            pout() << "Writing checkpoint file...\n";
            writeCheckpointFile();
        }
    }
    return steadyStateStop;
}
