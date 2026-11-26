"""Run all experiments and verify the saved tables using the bundled extracted data."""
import analysis
import early_sample
import turnaround
import manuscript_figures
import diagnostics
import verify
if __name__=="__main__":
    analysis.main()
    early_sample.main()
    turnaround.main()
    diagnostics.main()
    manuscript_figures.main()
    verify.main()
