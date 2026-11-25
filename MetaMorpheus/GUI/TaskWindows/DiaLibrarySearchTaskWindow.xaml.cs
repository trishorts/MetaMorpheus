using EngineLayer;
using GuiFunctions;
using GuiFunctions.Util;
using MassSpectrometry;
using Nett;
using Omics.Digestion;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.ComponentModel;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Input;
using TaskLayer;

namespace MetaMorpheusGUI
{
    /// <summary>
    /// The window that adds or edits a <see cref="DiaLibrarySearchTask"/>. It shows only what the DIA search reads: its own
    /// parameters (through <see cref="TaskWindowDiaLibrarySearch"/>), the task name, threads and q-value threshold, and the
    /// digestion and modifications used for protein groups and for a predicted library.
    /// </summary>
    public partial class DiaLibrarySearchTaskWindow : Window
    {
        private readonly ObservableCollection<ModTypeForTreeViewModel> FixedModTypeForTreeViewObservableCollection = new();
        private readonly ObservableCollection<ModTypeForTreeViewModel> VariableModTypeForTreeViewObservableCollection = new();

        public DiaLibrarySearchTaskWindow(DiaLibrarySearchTask task)
        {
            InitializeComponent();
            TheTask = task ?? new DiaLibrarySearchTask(); // null when there is no saved default
            PopulateChoices();
            UpdateFieldsFromTask(TheTask);
            if (task == null)
            {
                SaveButton.Content = "Add the DIA Library Search Task";
            }
            SearchModifications.Timer.Tick += TextChangeTimerHandler;
            Closing += OnClosing;
        }

        internal DiaLibrarySearchTask TheTask { get; private set; }

        private void PopulateChoices()
        {
            // Custom needs the custom-fragmentation window, which a prediction model cannot use
            foreach (string dissociationType in GlobalVariables.AllSupportedDissociationTypes.Keys.Where(k => k != DissociationType.Custom.ToString()))
            {
                DissociationTypeComboBox.Items.Add(dissociationType);
            }
            foreach (DiaLibraryDecoys decoys in Enum.GetValues<DiaLibraryDecoys>())
            {
                LibraryDecoysComboBox.Items.Add(decoys);
            }
            foreach (Protease protease in ProteaseDictionary.Dictionary.Values)
            {
                ProteaseComboBox.Items.Add(protease);
            }
            foreach (string behavior in Enum.GetNames(typeof(InitiatorMethionineBehavior)))
            {
                InitiatorMethionineBehaviorComboBox.Items.Add(behavior);
            }
            foreach (var modsOfAType in GlobalVariables.AllModsKnown.GroupBy(b => b.ModificationType))
            {
                var fixedType = new ModTypeForTreeViewModel(modsOfAType.Key, false);
                var variableType = new ModTypeForTreeViewModel(modsOfAType.Key, false);
                FixedModTypeForTreeViewObservableCollection.Add(fixedType);
                VariableModTypeForTreeViewObservableCollection.Add(variableType);
                foreach (var mod in modsOfAType)
                {
                    fixedType.Children.Add(new ModForTreeViewModel(mod.ToString(), false, mod.IdWithMotif, false, fixedType));
                    variableType.Children.Add(new ModForTreeViewModel(mod.ToString(), false, mod.IdWithMotif, false, variableType));
                }
            }
            FixedModsTreeView.DataContext = FixedModTypeForTreeViewObservableCollection;
            VariableModsTreeView.DataContext = VariableModTypeForTreeViewObservableCollection;
        }

        private void UpdateFieldsFromTask(DiaLibrarySearchTask task)
        {
            CommonParameters common = task.CommonParameters;
            var shown = TaskWindowDiaLibrarySearch.Load(task.Parameters);
            FragmentToleranceTextBox.Text = shown.FragmentTolerancePpm;
            Ms1ToleranceTextBox.Text = shown.Ms1TolerancePpm;
            ClassifierSeedTextBox.Text = shown.ClassifierSeed;
            LibraryChargesTextBox.Text = shown.LibraryPrecursorCharges;
            LibraryNceTextBox.Text = shown.LibraryNce;
            LibraryDecoysComboBox.SelectedItem = shown.LibraryDecoys;

            OutputFileNameTextBox.Text = common.TaskDescriptor;
            MaxThreadsTextBox.Text = common.MaxThreadsToUsePerFile.ToString(CultureInfo.InvariantCulture);
            FilterThresholdTextBox.Text = common.FilterThreshold.ToString(CultureInfo.InvariantCulture);
            DissociationTypeComboBox.SelectedItem = common.DissociationType.ToString();

            if (common.DigestionParams is DigestionParams digestion)
            {
                ProteaseComboBox.SelectedItem = TaskWindowSearchMode.ProteaseToShow(digestion);
                string searchModeWarning = TaskWindowSearchMode.PreserveWarning(digestion);
                SearchModeWarningTextBlock.Text = searchModeWarning;
                SearchModeWarningTextBlock.Visibility = searchModeWarning == null ? Visibility.Collapsed : Visibility.Visible;
                InitiatorMethionineBehaviorComboBox.SelectedIndex = (int)digestion.InitiatorMethionineBehavior;
            }
            IDigestionParams loaded = common.DigestionParams;
            MissedCleavagesTextBox.Text = loaded.MaxMissedCleavages == int.MaxValue ? "" : loaded.MaxMissedCleavages.ToString(CultureInfo.InvariantCulture);
            MinPeptideLengthTextBox.Text = loaded.MinLength.ToString(CultureInfo.InvariantCulture);
            MaxPeptideLengthTextBox.Text = loaded.MaxLength == int.MaxValue ? "" : loaded.MaxLength.ToString(CultureInfo.InvariantCulture);
            MaxModsPerPeptideTextBox.Text = loaded.MaxMods.ToString(CultureInfo.InvariantCulture);

            ShowMods(common.ListOfModsFixed, FixedModTypeForTreeViewObservableCollection);
            ShowMods(common.ListOfModsVariable, VariableModTypeForTreeViewObservableCollection);
        }

        /// <summary>Ticks each of the task's modifications in a tree, adding any the GUI does not know, as the other task windows do.</summary>
        private static void ShowMods(IEnumerable<(string, string)> mods, ObservableCollection<ModTypeForTreeViewModel> tree)
        {
            foreach (var (type, name) in mods)
            {
                var theModType = tree.FirstOrDefault(b => b.DisplayName.Equals(type));
                if (theModType == null)
                {
                    theModType = new ModTypeForTreeViewModel(type, true);
                    tree.Add(theModType);
                }
                var theMod = theModType.Children.FirstOrDefault(b => b.ModName.Equals(name));
                if (theMod != null)
                    theMod.Use = true;
                else
                    theModType.Children.Add(new ModForTreeViewModel("UNKNOWN MODIFICATION!", true, name, true, theModType));
            }
            foreach (var modType in tree)
            {
                modType.VerifyCheckState();
            }
        }

        private static List<(string, string)> TickedMods(ObservableCollection<ModTypeForTreeViewModel> tree) =>
            tree.SelectMany(type => type.Children.Where(b => b.Use).Select(b => (b.Parent.DisplayName, b.ModName))).ToList();

        private void SaveButton_Click(object sender, RoutedEventArgs e)
        {
            const string fieldNotUsed = "1";
            if (!TaskValidator.CheckTaskSettingsValidity(fieldNotUsed, fieldNotUsed, null, MissedCleavagesTextBox.Text,
                    fieldNotUsed, MinPeptideLengthTextBox.Text, MaxPeptideLengthTextBox.Text, MaxThreadsTextBox.Text, fieldNotUsed,
                    fieldNotUsed, fieldNotUsed, fieldNotUsed, fieldNotUsed, fieldNotUsed, fieldNotUsed, null, null, fieldNotUsed,
                    MaxModsPerPeptideTextBox.Text, fieldNotUsed, FilterThresholdTextBox.Text, null))
            {
                return;
            }
            var shown = new TaskWindowDiaLibrarySearch.Fields(FragmentToleranceTextBox.Text, Ms1ToleranceTextBox.Text,
                ClassifierSeedTextBox.Text, LibraryChargesTextBox.Text, LibraryNceTextBox.Text, (DiaLibraryDecoys)LibraryDecoysComboBox.SelectedItem);
            if (TaskWindowDiaLibrarySearch.Error(shown) is string error)
            {
                MessageBox.Show(error, "DIA Library Search Task", MessageBoxButton.OK, MessageBoxImage.Warning);
                return;
            }
            var listOfModsVariable = TickedMods(VariableModTypeForTreeViewObservableCollection);
            if (!TaskValidator.VariableModCheck(listOfModsVariable))
            {
                return;
            }

            var protease = (DigestionAgent)ProteaseComboBox.SelectedItem;
            int maxMissedCleavages = string.IsNullOrEmpty(MissedCleavagesTextBox.Text) ? int.MaxValue : int.Parse(MissedCleavagesTextBox.Text, CultureInfo.InvariantCulture);
            int minPeptideLength = int.Parse(MinPeptideLengthTextBox.Text, NumberStyles.Any, CultureInfo.InvariantCulture);
            int maxPeptideLength = string.IsNullOrEmpty(MaxPeptideLengthTextBox.Text) ? int.MaxValue : int.Parse(MaxPeptideLengthTextBox.Text, NumberStyles.Any, CultureInfo.InvariantCulture);
            int maxModsPerPeptide = int.Parse(MaxModsPerPeptideTextBox.Text, CultureInfo.InvariantCulture);
            int threads = int.Parse(MaxThreadsTextBox.Text, CultureInfo.InvariantCulture);
            CommonParameters loaded = TheTask.CommonParameters;

            // pass only what this window shows; the search mode and every other digestion setting are kept from the loaded task
            var digestion = TaskWindowDigestionParams.ForSave(loaded.DigestionParams, protease.Name, maxMissedCleavages, minPeptideLength,
                maxPeptideLength, loaded.DigestionParams.MaxModificationIsoforms, maxModsForPeptides: maxModsPerPeptide,
                initiatorMethionineBehavior: InitiatorMethionineBehaviorComboBox.SelectedIndex >= 0
                    ? (InitiatorMethionineBehavior)InitiatorMethionineBehaviorComboBox.SelectedIndex
                    : null);

            TheTask.Parameters = TaskWindowDiaLibrarySearch.ForSave(shown);
            TheTask.CommonParameters = new CommonParameters(
                taskDescriptor: OutputFileNameTextBox.Text != "" ? OutputFileNameTextBox.Text : "DiaLibrarySearchTask",
                maxThreadsToUsePerFile: threads > 0 && threads <= Environment.ProcessorCount ? threads : new CommonParameters().MaxThreadsToUsePerFile,
                digestionParams: digestion,
                dissociationType: GlobalVariables.AllSupportedDissociationTypes[DissociationTypeComboBox.SelectedItem.ToString()],
                filterThreshold: double.Parse(FilterThresholdTextBox.Text, NumberStyles.Float, CultureInfo.InvariantCulture),
                listOfModsFixed: TickedMods(FixedModTypeForTreeViewObservableCollection),
                listOfModsVariable: listOfModsVariable,
                // read by database loading, which a protein database goes through; the window has no box for them
                addTruncations: loaded.AddTruncations,
                minVariantDepth: loaded.MinVariantDepth,
                maxHeterozygousVariants: loaded.MaxHeterozygousVariants);
            DialogResult = true;
        }

        private void CancelButton_Click(object sender, RoutedEventArgs e)
        {
            DialogResult = false;
        }

        private void KeyPressed(object sender, KeyEventArgs e)
        {
            if (e.Key == Key.Return)
            {
                SaveButton_Click(sender, e);
            }
            else if (e.Key == Key.Escape)
            {
                CancelButton_Click(sender, e);
            }
        }

        private void SaveAsDefault_Click(object sender, RoutedEventArgs e)
        {
            SaveButton_Click(sender, e);
            Toml.WriteFile(TheTask, Path.Combine(GlobalVariables.DataDir, "DefaultParameters", "DiaLibrarySearchTaskDefault.toml"), MetaMorpheusTask.tomlConfig);
        }

        private void TextChanged_Fixed(object sender, TextChangedEventArgs args)
        {
            SearchModifications.SetTimer();
            SearchModifications.FixedSearch = true;
        }

        private void TextChanged_Var(object sender, TextChangedEventArgs args)
        {
            SearchModifications.SetTimer();
            SearchModifications.VariableSearch = true;
        }

        private void TextChangeTimerHandler(object sender, EventArgs e)
        {
            if (SearchModifications.FixedSearch)
            {
                SearchModifications.FilterTree(SearchFixMod, FixedModsTreeView, FixedModTypeForTreeViewObservableCollection);
                SearchModifications.FixedSearch = false;
            }
            if (SearchModifications.VariableSearch)
            {
                SearchModifications.FilterTree(SearchVarMod, VariableModsTreeView, VariableModTypeForTreeViewObservableCollection);
                SearchModifications.VariableSearch = false;
            }
        }

        private void OnClosing(object sender, CancelEventArgs e)
        {
            // a closed window must not stay on the shared timer, or the next tick finds it gone
            SearchModifications.Timer.Tick -= TextChangeTimerHandler;
        }
    }
}
