package cli

import (
	"fmt"

	"github.com/compgenlab/gtstore/internal/buildinfo"
	"github.com/spf13/cobra"
)

// Set via -ldflags at build time.
var (
	Version = "dev"
	GitHash = ""
)

func versionString() string {
	if GitHash != "" {
		return fmt.Sprintf("%s (%s)", Version, GitHash)
	}
	return Version
}

var rootCmd = &cobra.Command{
	Use:   "gtstore",
	Short: "Sparse genotype store for cohort VCFs and gVCFs",
	Long: `gtstore converts cohort VCFs into a sparse, queryable genotype store and
answers genotype queries against a store or a VCF alike.

  gtstore convert        --out cohort/chr1 chr1.vcf.gz      VCF -> volume
  gtstore create         cohort/                            volumes -> gtstore
  gtstore refblock-calc  --roster cohort --out frags/S1 S1.g.vcf.gz
  gtstore refblock-attach cohort/ frags/                    gVCF blocks -> store
  gtstore query          cohort --variant chr1:1000:A:T     genotypes
  gtstore summary        cohort                             what a store holds`,
	// Silence usage only once we reach a command's RunE. Cobra validates args
	// and flags *before* PersistentPreRun, so parsing errors (e.g. a missing
	// argument) still print usage, while errors during execution do not.
	PersistentPreRun: func(cmd *cobra.Command, args []string) {
		cmd.SilenceUsage = true
	},
}

// Run executes the root command and reports whether it failed, leaving the exit
// to the caller so that deferred cleanup -- notably stopping a CPU profile and
// closing its file -- still runs.
func Run() error { return rootCmd.Execute() }

func init() {
	buildinfo.Version = Version
	buildinfo.GitHash = GitHash

	rootCmd.Version = versionString()
	rootCmd.CompletionOptions.HiddenDefaultCmd = true

	// Cobra's default help, plus a "Since:" line from each command's "since"
	// annotation, plus a version footer.
	const helpBody = `{{with (or .Long .Short)}}{{. | trimTrailingWhitespaces}}

{{end}}{{if .Annotations}}{{with (index .Annotations "since")}}Since: {{.}}

{{end}}{{end}}{{if or .Runnable .HasSubCommands}}{{.UsageString}}{{end}}`
	rootCmd.SetHelpTemplate(helpBody + "\ngtstore " + versionString() +
		" https://github.com/compgenlab/gtstore\n")

	InitCmd(rootCmd)
}
