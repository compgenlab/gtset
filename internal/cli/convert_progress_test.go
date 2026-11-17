package cli

import (
	"bytes"
	"fmt"
	"os"
	"path/filepath"
	"strings"
	"testing"
	"time"

	"github.com/compgenlab/cghts/gtstore"
)

// progressFixture writes a VCF of n records on chr1, step bases apart,
// declaring chr1's length when length > 0.
func progressFixture(t *testing.T, n, step int, length int64) string {
	t.Helper()
	var b strings.Builder
	b.WriteString("##fileformat=VCFv4.2\n")
	if length > 0 {
		fmt.Fprintf(&b, "##contig=<ID=chr1,length=%d>\n", length)
	} else {
		b.WriteString("##contig=<ID=chr1>\n")
	}
	b.WriteString("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">\n")
	b.WriteString("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"DP\">\n")
	b.WriteString("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n")
	for i := 0; i < n; i++ {
		fmt.Fprintf(&b, "chr1\t%d\t.\tA\tG\t.\tPASS\t.\tGT:DP\t0/1:30\t0/0:30\n", 1+i*step)
	}
	p := filepath.Join(t.TempDir(), "in.vcf")
	if err := os.WriteFile(p, []byte(b.String()), 0o644); err != nil {
		t.Fatal(err)
	}
	return p
}

// fastClock makes every reading of the clock 301 seconds later than the last,
// so the five-minute interval has always passed when the converter looks.
func fastClock(t *testing.T) {
	t.Helper()
	at := time.Date(2026, 1, 1, 0, 0, 0, 0, time.UTC)
	old := progressClock
	progressClock = func() time.Time { at = at.Add(301 * time.Second); return at }
	t.Cleanup(func() { progressClock = old })
}

func progressLines(out string) []string {
	var lines []string
	for _, l := range strings.Split(out, "\n") {
		if strings.Contains(l, "records ") && strings.Contains(l, "/s)") {
			lines = append(lines, l)
		}
	}
	return lines
}

func TestConvertReportsProgressOnAClock(t *testing.T) {
	fastClock(t)
	in := progressFixture(t, 300, 30, 10_000) // 64 records move 1,920 bases: >1% of the contig
	out := runVcf(t, "convert", "-v", "--out", filepath.Join(t.TempDir(), "s"), in)
	lines := progressLines(out)
	if len(lines) == 0 {
		t.Fatalf("no progress lines in:\n%s", out)
	}
	for _, want := range []string{"chr1:", "/ 10,000 (", "%)", "eta ~", "calls "} {
		if !strings.Contains(lines[0], want) {
			t.Errorf("progress line lacks %q: %s", want, lines[0])
		}
	}
	if !strings.Contains(out, "chr1: starting at 1 of 10,000") {
		t.Errorf("no chromosome start line in:\n%s", out)
	}
	if !strings.Contains(out, "records/s") || !strings.Contains(out, "genotypes/s") {
		t.Errorf("summary lacks throughput:\n%s", out)
	}
}

// Five minutes that moved two bases is not news.
func TestConvertProgressNeedsRealMovement(t *testing.T) {
	fastClock(t)
	in := progressFixture(t, 300, 2, 1_000_000_000) // under 1,000 records and far under 1% of the contig
	out := runVcf(t, "convert", "-v", "--out", filepath.Join(t.TempDir(), "s"), in)
	if lines := progressLines(out); len(lines) != 0 {
		t.Errorf("progress printed without real movement: %q", lines)
	}
}

// Without a declared length there is nothing to take a percentage of or
// estimate from; the line says what it can and the record threshold still holds.
func TestConvertProgressWithoutContigLength(t *testing.T) {
	fastClock(t)
	in := progressFixture(t, 1_300, 2, 0)
	out := runVcf(t, "convert", "-v", "--out", filepath.Join(t.TempDir(), "s"), in)
	lines := progressLines(out)
	if len(lines) != 1 {
		t.Fatalf("want exactly one line (after 1,000 records), got %q", lines)
	}
	if strings.Contains(lines[0], "%") || strings.Contains(lines[0], "eta") {
		t.Errorf("line claims a percentage or estimate it cannot know: %s", lines[0])
	}
}

func TestConvertProgressSilentWithoutVerbose(t *testing.T) {
	fastClock(t)
	in := progressFixture(t, 300, 30, 10_000)
	out := runVcf(t, "convert", "--out", filepath.Join(t.TempDir(), "s"), in)
	if lines := progressLines(out); len(lines) != 0 {
		t.Errorf("progress without -v: %q", lines)
	}
}

type flushCounter struct {
	bytes.Buffer
	flushes int
}

func (f *flushCounter) Flush() error { f.flushes++; return nil }

// A cluster job's log is only useful if a line is in it while the job runs.
func TestProgressLinesAreFlushed(t *testing.T) {
	w, err := gtstore.NewWriter(filepath.Join(t.TempDir(), "s"), gtstore.WriterOpts{Samples: []string{"S1"}, MinDP: 10})
	if err != nil {
		t.Fatal(err)
	}
	defer w.Discard()
	out := &flushCounter{}
	at := time.Date(2026, 1, 1, 0, 0, 0, 0, time.UTC)
	c := &parquetConverter{w: w, verbose: true, progress: out,
		now: func() time.Time { at = at.Add(301 * time.Second); return at }}
	c.curChrom = "chr1"
	c.note("chr1", 1)
	if out.flushes != 1 {
		t.Fatalf("start line: %d flushes, want 1", out.flushes)
	}
	c.nRecords = 64 * 20 // past the record threshold, on a clock check
	c.tick(5_000)
	if out.flushes != 2 || !strings.Contains(out.String(), "chr1:5,000") {
		t.Errorf("progress line: %d flushes, output %q", out.flushes, out.String())
	}
}

func TestRoughDuration(t *testing.T) {
	for d, want := range map[time.Duration]string{
		30 * time.Second:                "<1m",
		47*time.Minute + 40*time.Second: "48m",
		time.Hour + 12*time.Minute:      "1h10m",
		2*time.Hour + 20*time.Minute:    "2h15m",
		26*time.Hour + 2*time.Minute:    "26h0m",
	} {
		if got := roughDuration(d); got != want {
			t.Errorf("roughDuration(%s) = %q, want %q", d, got, want)
		}
	}
}
