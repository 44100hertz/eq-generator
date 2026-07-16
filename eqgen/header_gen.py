"""
Shared C header generation for eq_coeffs.h.

Used by both the wire CLI (live PipeWire filter) and the export CLI
(ESP32 firmware export).  Each caller supplies its own define lines
so it can control exactly which EQGEN_* constants are emitted and
in what order — the shared function handles the common formatting:
header comment, coefficient arrays, and the eqgen_get_coeffs /
eqgen_get_fs inline helpers.
"""


def _format_coeff_array(name, biquads, bands):
    """Format one coefficient array as a list of C source lines.

    Each biquad writes 5 float coefficients with a per-band comment.
    """
    lines = [f"static const float {name}[{len(bands) * 5}] = {{"]
    for i, bc in enumerate(biquads):
        band = bands[i]
        lines.append(
            f"    {bc.b0:>15.9f}f, {bc.b1:>15.9f}f, {bc.b2:>15.9f}f, "
            f"{bc.a1:>15.9f}f, {bc.a2:>15.9f}f,"
            f"  // [{i}] f0={band['f0']:6.1f} Hz gain={band['gain_db']:+5.1f} dB"
        )
    lines.append("};")
    lines.append("")
    return lines


def generate_c_header(
    biquads_44, bands_44,
    biquads_48, bands_48,
    *,
    header_lines,
    define_lines,
) -> str:
    """Build the full eq_coeffs.h content.

    Parameters
    ----------
    biquads_44, bands_44 : 44.1 kHz biquad coefficients and band metadata.
    biquads_48, bands_48 : 48 kHz biquad coefficients and band metadata.
    header_lines : list[str]
        File-top comment lines (without ``//`` prefix).
    define_lines : list[str]
        Raw lines for the define block — emitted verbatim between
        ``#pragma once`` and the coefficient arrays.  Each caller
        assembles exactly the ``EQGEN_*`` defines it needs.

    Returns
    -------
    str
        Complete ``eq_coeffs.h`` file content.
    """
    lines = []

    # ── Header comment ────────────────────────────────────────────
    for hl in header_lines:
        lines.append(f"// {hl}")
    lines.append("")
    lines.append("#pragma once")
    lines.append("")

    # ── Caller-supplied defines ───────────────────────────────────
    for dl in define_lines:
        lines.append(dl)

    # ── Coefficient arrays ────────────────────────────────────────
    if len(bands_44) > 0:
        lines.extend(_format_coeff_array("eqgen_coeffs_44100", biquads_44, bands_44))
        lines.extend(_format_coeff_array("eqgen_coeffs_48000", biquads_48, bands_48))

    # ── Runtime selectors ─────────────────────────────────────────
    lines.append("/** Select coefficient array for the given sample rate. */")
    lines.append("static inline const float *eqgen_get_coeffs(int sample_rate) {")
    lines.append("    return (sample_rate == 48000) ? eqgen_coeffs_48000 : eqgen_coeffs_44100;")
    lines.append("}")
    lines.append("")
    lines.append("/** Select nominal Fs for the given sample rate. */")
    lines.append("static inline int eqgen_get_fs(int sample_rate) {")
    lines.append("    return (sample_rate == 48000) ? EQGEN_FS_48000 : EQGEN_FS_44100;")
    lines.append("}")

    return "\n".join(lines) + "\n"
