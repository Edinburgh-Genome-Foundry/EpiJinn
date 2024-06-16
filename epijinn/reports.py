# Copyright 2024 Edinburgh Genome Foundry, University of Edinburgh
#
# This file is part of EpiJinn.
#
# EpiJinn is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#
# EpiJinn is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with Ediacara. If not, see <https://www.gnu.org/licenses/>.

from datetime import datetime
import os

from pdf_reports import (
    add_css_class,
    dataframe_to_html,
    pug_to_html,
    style_table_rows,
    write_report,
)
import pdf_reports.tools as pdf_tools


from .version import __version__

THIS_PATH = os.path.dirname(os.path.realpath(__file__))
ASSETS_PATH = os.path.join(THIS_PATH, "report_assets")
BEDMETHYLITEMGROUP_REPORT_TEMPLATE = os.path.join(ASSETS_PATH, "report.pug")
STYLESHEET = os.path.join(ASSETS_PATH, "report_style.css")


def end_pug_to_html(template, **context):
    now = datetime.now().strftime("%Y-%m-%d")
    defaults = {
        "sidebar_text": "Generated on %s by EpiJinn (version %s)" % (now, __version__),
        "epijinn_logo_url": os.path.join(ASSETS_PATH, "imgs", "epijinn.png"),
    }
    for k in defaults:
        if k not in context:
            context[k] = defaults[k]
    return pug_to_html(template, **context)


def write_bedmethylitemgroup_report(target, bedmethylitemgroup):
    """Write a methylation analysis report with a PDF summary.


    **Parameters**

    **target**
    > Path for PDF file.

    **bedmethylitemgroup**
    > `BedmethylItemGroup` instance.
    """
    for bedmethylitem in bedmethylitemgroup.bedmethylitems:
        bedmethylitem.bedmethylitem_figure_data = pdf_tools.figure_data(
            bedmethylitem.fig, fmt="svg"
        )
        for bedresult in bedmethylitem.results:
            bedresult.bed_pdf = dataframe_to_html(
                bedresult.bed, extra_classes=("definition",)
            )
            bedresult.figure_data = pdf_tools.figure_data(bedresult.plot, fmt="svg")

    html = end_pug_to_html(
        BEDMETHYLITEMGROUP_REPORT_TEMPLATE, bedmethylitemgroup=bedmethylitemgroup
    )
    write_report(html, target, extra_stylesheets=(STYLESHEET,))
