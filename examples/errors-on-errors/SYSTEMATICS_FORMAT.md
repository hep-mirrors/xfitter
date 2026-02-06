# Errors-on-Errors (EoE) Systematics Configuration

This example demonstrates the **embedded epsilon format** for configuring Errors-on-Errors (EoE) corrections in xFitter.

## New Format (Embedded Epsilon)

The new format embeds epsilon values directly in the source specification string using the `@eps=` suffix:

```fortran
&Systematics
  ListOfSources =
    'sysHZComb1036:N@eps=0.7',
    'sysHZComb1042:N@eps=0.7',
    'proc_tb21:E@eps=0.5',
    'sysHZComb1041:N@eps=0.3'

  n_iterations = 2
  Enable_Bartlett = .true.
&End
```

## Syntax

Each source string follows the pattern:
```
'sourcename:modifiers@eps=value'
```

Where:
- **sourcename**: Name of the systematic source (must match data file sources)
- **:modifiers**: Optional form and scaling modifiers (see below)
- **@eps=value**: Epsilon value for EoE correction (must be >= 0)

## Supported Modifiers

### Chi2 Form Modifiers (EoE compatible: N, E only)

- **:N** - Nuisance parameters (Hessian method) - **EoE compatible**
- **:E** - External (minimizer method) - **EoE compatible**
- **:C** - Covariance matrix method - **NOT compatible with EoE**
- **:O** - Offset method - **NOT compatible with EoE**

### Scaling Modifiers

- **:M** - Multiplicative/Linear scaling
- **:A** - Additive scaling (no rescaling)
- **:P** - Poisson scaling

### Type Modifiers

- **:D** - Data systematic (default)
- **:T** - Theory systematic

## Examples

### Basic EoE Configuration

```fortran
ListOfSources =
  'sysLumi:N@eps=0.7',
  'sysJES:N@eps=0.5'
```

### Multiple Modifiers with EoE

```fortran
ListOfSources =
  'sysLumi:N:M@eps=0.7',           ! Nuisance, multiplicative, eps=0.7
  'sysTheory:T:N@eps=0.5',         ! Theory, nuisance, eps=0.5
  'sysBeam:E:A@eps=0.3'            ! External, additive, eps=0.3
```

## Important Notes

1. **Form Compatibility**: EoE only works with `:N` (Nuisance) and `:E` (External) forms. Using `@eps=` with `:C` or `:O` will trigger a warning and the epsilon value will be ignored.

2. **Epsilon Requirements**:
   - Must be >= 0
   - No parallel array alignment issues - each source is self-contained

3. **Modifiers Before @eps**: Always place form/scaling modifiers (`:N`, `:M`, etc.) before the `@eps=` suffix.

## Additional EoE Parameters

Beyond the epsilon values, you can configure:

- **n_iterations**: Number of EoE iterations (default: 2)
- **Enable_Bartlett**: Enable Bartlett corrections (default: .true.)

```fortran
&Systematics
  ListOfSources = 'sysA:N@eps=0.7', 'sysB:N@eps=0.5'
  n_iterations = 2
  Enable_Bartlett = .true.
&End
```

## References

For more details on Errors-on-Errors methodology, see:
- Errors-on-errors: Eur.Phys.J.C 85 (2025)
- Higher-order Asymtotic corrections: Eur.Phys.J.C 83 (2023)
