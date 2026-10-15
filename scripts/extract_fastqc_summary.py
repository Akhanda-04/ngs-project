
from pathlib import Path
from html.parser import HTMLParser

REPORT_DIRS = [
    Path("qc_reports/pre_trim"),
    Path("qc_reports/post_trim"),
]


class FastQCParser(HTMLParser):
    def __init__(self):
        super().__init__()
        self.in_table = False
        self.in_cell = False
        self.current_row = []
        self.rows = []
        self.current_text = []

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if tag == "table":
            self.in_table = True
        elif self.in_table and tag == "tr":
            self.current_row = []
        elif self.in_table and tag in ("td", "th"):
            self.in_cell = True
            self.current_text = []

    def handle_data(self, data):
        if self.in_table and self.in_cell:
            self.current_text.append(data.strip())

    def handle_endtag(self, tag):
        if self.in_table and tag in ("td", "th") and self.in_cell:
            value = " ".join(x for x in self.current_text if x)
            self.current_row.append(value)
            self.in_cell = False
        elif self.in_table and tag == "tr":
            if self.current_row:
                self.rows.append(self.current_row)
        elif tag == "table":
            self.in_table = False


for directory in REPORT_DIRS:
    for report in sorted(directory.glob("*_fastqc.html")):
        parser = FastQCParser()
        parser.feed(report.read_text(errors="replace"))

        print(f"\n{'=' * 65}")
        print(f"REPORT: {report}")
        print("=" * 65)

        for row in parser.rows:
            if len(row) >= 2:
                print(f"{row[0]}: {row[1]}")
